//
// $Id$
//
//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> manchester.ac.uk>
//
// Copyright (C) 2013  CADET Laboratory for Medical Bioinformatics, University of Manchester, UK
//
// This file is part of seaMass.
//
// seaMass is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// seaMass is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with seaMass.  If not, see <http://www.gnu.org/licenses/>.
//

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <limits>
#include <math.h>

#include <omp.h>
#include <SpatialIndex.h>


using namespace std;
using namespace SpatialIndex;


struct MyPair
{
	id_type id;
	Region* r;

	MyPair() {}

	MyPair(const id_type& _id, Region* _r) :
		id(_id),
		r(_r) {}

	void operator=(const MyPair& pair)
	{
		id = pair.id;
		r = pair.r;
	}
};

// yes, this is reversed
inline bool operator<(const MyPair& a, const MyPair& b)
{
	return a.r->getLow(2) > b.r->getLow(2);
}

// if all the nodes were internally sorted by intensity
// we could use a quicker sorting procedure than priority_queue
class MyQueryStrategy : public SpatialIndex::IQueryStrategy
{
private:
	const IShape& query;
	string basename;

	priority_queue<MyPair> node_q;
	priority_queue<MyPair> data_q;
	vector<MyPair> chunk;
	int chunk_index;
	int index;

	double start;

public:
	MyQueryStrategy(const IShape& _query,
		            string& _basename,
					int _chunk_size = 4096) : 
	  query(_query),
	  basename(_basename),
	  chunk(_chunk_size),
	  chunk_index(0),
	  index(0)
	{
		ostringstream oss; oss << basename << "_" << setfill('0') << setw(5) << chunk_index++ << "_0.out";
		ofstream ofs(oss.str().c_str());
		start = omp_get_wtime();
	}

	void getNextEntry(const IEntry& entry, id_type& nextEntry, bool& hasNext)
	{
		if (!&entry)
		{
			hasNext = false;
			return;
		}

		const INode* n = dynamic_cast<const INode*>(&entry);

		IShape* ps;
		n->getShape(&ps);
		Region* pr = dynamic_cast<Region*>(ps);

		for (uint32_t cChild = 0; cChild < n->getChildrenCount(); cChild++)
		{
			IShape* ps;
			n->getChildShape(cChild, &ps);
			Region* r = dynamic_cast<Region*>(ps);
			if (query.intersectsShape(*r))
			{
				if (n->getLevel() == 0)
					data_q.push(MyPair(n->getChildIdentifier(cChild), r));
				else
					node_q.push(MyPair(n->getChildIdentifier(cChild), r));
			}
		}

		while (node_q.empty() && !data_q.empty() ||
			   !data_q.empty() && data_q.top().r->getLow(2) <= node_q.top().r->getLow(2))
		{
			visitData(data_q.top());
			data_q.pop();
		}

		if (!node_q.empty())
		{
			Region* r = node_q.top().r;
			nextEntry = node_q.top().id; node_q.pop();
			hasNext = true;
			delete r;
		}
		else
		{
			hasNext = false;
		}
	}

	void visitData(const MyPair& pair)
	{
		chunk[index] = pair;

		if (++index == chunk.size())		
		{
			ostringstream oss; oss << basename << "_" << setfill('0') << setw(5) << chunk_index << "_" << omp_get_wtime() - start << ".out";
			ofstream ofs(oss.str().c_str());
			for (size_t i = 0; i < chunk.size(); i++)
			{
				double w = 0.25 * (chunk[i].r->getHigh(0) - chunk[i].r->getLow(0));
				double h = 0.25 * (chunk[i].r->getHigh(1) - chunk[i].r->getLow(1));
				int x = chunk[i].r->getLow(0) / w + 3;
				int y = chunk[i].r->getLow(1) / h + 3;
				int lx = log(w) / log(2.0);
				int ly = log(h) / log(2.0);
				double intensity = -chunk[i].r->getLow(2);
				
				ofs << intensity << ":" << lx << "," << ly << ":" << x << "," << y << endl;
			
				delete chunk[i].r;
			}
			cout << "." << flush;

			index = 0;
			chunk_index++;
		}
	}
};


int main(int argc, char *argv[])
{
    cout << endl;
    cout << "seaMass Viz Client - Copyright (C) 2014 - biospI Laboratory" << endl;
    cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
    cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;  
	cout << endl;	
	if (argc != 8)
	{
		cout << "Usage" << endl;
		cout << "-----" << endl;
		cout << "viz_client <in_idx> <mz_min> <mz_max> <rt_min> <rt_max> <out_w> <out_h>" << endl;
		cout << endl;
		cout << "<in_dir>:  Input idx file" << endl;
		cout << "<mz_min>:  Minimum m/z to display" << endl;
		cout << "<mz_max>:  Maxmimum m/z to display" << endl;
		cout << "<rt_min>:  Minimum retention time to display" << endl;
		cout << "<rt_max>:  Maximum retention time to display" << endl;
		cout << "<out_w>:   Output width in pixels" << endl;
		cout << "<out_h>:   Output height in pixels" << endl;
		return 0;
	}
	string in_dir(argv[1]);
	double mz_min = atof(argv[2]);
	double mz_max = atof(argv[3]);
	double rt_min = atof(argv[4]);
	double rt_max = atof(argv[5]);
	int out_w  = atoi(argv[6]);
	int out_h  = atoi(argv[7]);

	int lastdot = in_dir.find_last_of("."); 
	string basename = (lastdot == string::npos) ? in_dir : in_dir.substr(0, lastdot);

	IStorageManager* diskfile = StorageManager::loadDiskStorageManager(basename);
	// this will try to locate and open an already existing storage manager.

	StorageManager::IBuffer* file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);
	// applies a main memory random buffer on top of the persistent storage manager
	// (LRU buffer, etc can be created the same way).

	// If we need to open an existing tree stored in the storage manager, we only
	// have to specify the index identifier as follows
	ISpatialIndex* tree = RTree::loadRTree(*file, 1);

	double low[3] = {
		mz_min * 60/1.0033548378,
		rt_min,
		-numeric_limits<double>::max()
	};
	double high[3] = {
		mz_max * 60/1.0033548378,
		rt_max,
		0.0
	};
	Region r(low, high, 3);

	MyQueryStrategy qs(r, basename);
	tree->queryStrategy(qs);

    cout << endl << endl;

	return 0;
}
