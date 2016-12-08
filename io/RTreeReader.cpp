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

#include "RTreeReader.hpp"

#include <limits>
#include <queue>
#include <math.h>


using namespace SpatialIndex;
using namespace std;


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
	return a.r->getLow(a.r->m_dimension - 1) > b.r->getLow(a.r->m_dimension - 1);
}


// if all the nodes were internally sorted by intensity
// we could use a quicker sorting procedure than priority_queue
class MyQueryStrategy : public SpatialIndex::IQueryStrategy
{
private:
	const IShape& query;

	priority_queue<MyPair> node_q;
	priority_queue<MyPair> data_q;

	SeamassCore::Output& output;

public:
	MyQueryStrategy(const IShape& _query, SeamassCore::Output& _output) :
		query(_query),
		output(_output)
	{
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
			MyPair d = data_q.top();
			data_q.pop();
			visitData(d, !(node_q.empty() && data_q.size() == 1));
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

	void visitData(const MyPair& pair, bool hasNext)
	{
		if (hasNext)
		{
			double area = 1.0;
			for (ii i = 0; i < pair.r->m_dimension - 1; i++)
			{
				double width = 0.25 * (pair.r->getHigh(i) - pair.r->getLow(i));
				output.offsets[i].push_back((int)(pair.r->getLow(i) / width) + 3);
				output.scales[i].push_back((int)(-log(width) / log(2.0)));
				area *= width;
			}
			output.weights.push_back((fp)-pair.r->getLow(pair.r->m_dimension - 1) / area);
		}
	}
};


RTreeReader::
RTreeReader(const string& filename)
{
	// this will try to locate and open an already existing storage manager.
	string fn = filename;
	diskfile = StorageManager::loadDiskStorageManager(fn);

	// applies a main memory random buffer on top of the persistent storage manager
	// (LRU buffer, etc can be created the same way).
	file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);

	// If we need to open an existing tree stored in the storage manager, we only
	// have to specify the index identifier as follows
	tree = RTree::loadRTree(*file, 1);
}


RTreeReader::
~RTreeReader()
{
	delete tree;
	delete file;
	delete diskfile;
}


// only supports cm with dimension 2 at present
void
RTreeReader::
read(SeamassCore::Output& output)
{
	ii dimensions = output.baselineExtent.size();
	output.scales.resize(dimensions);
	output.offsets.resize(dimensions);


	// todo: for setting an actual region of interest
	/*double low[3] = {
		mz_min * 60 / 1.0033548378,
		rt_min,
		-numeric_limits<double>::max()
	};
	double high[3] = {
		mz_max * 60 / 1.0033548378,
		rt_max,
		0.0
	};*/

	//everything for now
	vector<double> low(dimensions + 1);
	vector<double> high(dimensions + 1);
	for (ii i = 0; i < dimensions; i++)
	{
		low[i] = 0.0;
		high[i] = numeric_limits<double>::max();
	}
	low[dimensions] = -numeric_limits<double>::max();
	high[dimensions] = 0.0;

	Region r(low.data(), high.data(), dimensions + 1);
	MyQueryStrategy qs(r, output);
	tree->queryStrategy(qs);
}