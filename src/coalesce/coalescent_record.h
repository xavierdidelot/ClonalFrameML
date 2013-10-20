/*  Copyright 2013 Daniel Wilson.
 *
 *  coalescent_record.h
 *  Part of the coalesce library.
 *
 *  The myutils library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  The myutils library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with the myutils library. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _RECORD_H_
#define _RECORD_H_

class mt_node
{
	static int number;
public:
	/*Fixed once*/
	int id;

	/*Recyclable*/
	bool in_use;
	double time;
	double edge_time;
	double last_update;			// in a structured coalescent, the last time edge_time was updated
	class mt_node *ancestor;	//ptr to ancestor
	class mt_node *descendant[2];//vec of ptrs to descendant

public:
	mt_node() {};
	mt_node& initialize(const int id_in)	
	{
		id=id_in;
		recycle();
		return *this;
	}
	mt_node& recycle()
	{
		in_use=false;
		time=0.0;
		edge_time=0.0;
		last_update=0.0;
		ancestor=NULL;
		descendant[0]=NULL;
		descendant[1]=NULL;
		return *this;
	}
};

class marginal_tree
{
	int id;
	int k;
public:
	/*Fixed once*/
	int size;
	//class Control *con;			//ptr to con
	int n;
	class mt_node *node;		//vec of mt_node's
	
	/*Recyclable*/
	int genotype;
	int next_free_node;
	int nco;
	
public:
	marginal_tree() {};
	/*marginal_tree& initialize(const int id_in, class Control *con_in)
	{
		id=id_in;
		con=con_in;
		size=con->nsamp+(con->nsamp-1);

		node=(class mt_node*) malloc((size_t) size*sizeof(class mt_node));
		int i;
		for(i=0;i<size;i++)
			node[i].initialize(i);	

		recycle();
		return *this;
	}*/
	marginal_tree& initialize(const int id_in, const int n_in)
	{
		id=id_in;
		n = n_in;
		size=n+(n-1);

		node=(class mt_node*) malloc((size_t) size*sizeof(class mt_node));
		int i;
		for(i=0;i<size;i++)
			node[i].initialize(i);	

		recycle();
		return *this;
	}
	marginal_tree& recycle()
	{
		int i;
		for(i=0;i<size;i++)
			node[i].recycle();
		genotype=-1;
		k=0;
		next_free_node=n;
		nco=0;
		return *this;
	}
	class mt_node* coalesce(double &time, int &id1, int &id2)
	{
		/**Additional error checking added 12/04/04						   **/
		if(id1==id2)
			error("coalesce(): cannot coalesce a node with itself");
		/**NB: if this error is triggered it is not because the same active**/
		/**nodes were chosen to coalesce with one another, but because they**/
		/**point to the same marginal node at that site, which is erroneous**/
		++nco;
		if(k==1)
		{
			error("coalesce(): called when no longer segregating");
		}
		node[next_free_node].in_use=true;
		node[next_free_node].time=time;
		node[id1].edge_time=time-node[id1].time;
		node[id2].edge_time=time-node[id2].time;
		node[id1].ancestor=&(node[next_free_node]);
		node[id2].ancestor=&(node[next_free_node]);
		node[next_free_node].descendant[0]=&(node[id1]);
		node[next_free_node].descendant[1]=&(node[id2]);

		//if(id==0)
		//	printf("Tree id's: New node: %d Coalescing nodes: %d and %d\n",node[next_free_node].id,node[id1].id,node[id2].id);

		++next_free_node;
		--k;
		return &node[next_free_node-1];
	}
	class mt_node* migrate_coalesce(double &time, int &id1, int &id2)
	{
		/**Additional error checking added 12/04/04						   **/
		if(id1==id2)
			error("coalesce(): cannot coalesce a node with itself");
		/**NB: if this error is triggered it is not because the same active**/
		/**nodes were chosen to coalesce with one another, but because they**/
		/**point to the same marginal node at that site, which is erroneous**/
		++nco;
		if(k==1)
		{
			error("coalesce(): called when no longer segregating");
		}
		node[next_free_node].in_use=true;
		node[next_free_node].time=time;
		node[next_free_node].edge_time=0.0;
		node[next_free_node].last_update=time;
		node[id1].ancestor=&(node[next_free_node]);
		node[id2].ancestor=&(node[next_free_node]);
		node[next_free_node].descendant[0]=&(node[id1]);
		node[next_free_node].descendant[1]=&(node[id2]);

		//if(id==0)
		//	printf("Tree id's: New node: %d Coalescing nodes: %d and %d\n",node[next_free_node].id,node[id1].id,node[id2].id);

		++next_free_node;
		--k;
		return &node[next_free_node-1];
	}
	class mt_node* add_base_node(double *time_in, int &id)
	{
		node[id].in_use=true;
		node[id].time=*time_in;
		k++;
		return &node[id];
	}
	int get_k() {return k;};
	void error(char* error_text)
	{
		printf("Error in marginal_tree::%s\n",error_text);
		printf("Exiting to system");
		exit(13);
	}

	double height() {
		return node[size-1].time;
	}
	double branch_length() {
		double result = 0.0;
		int i;
		for(i=0;i<size-1;i++) result += node[i].edge_time;
		return result;
	}
};

#endif
