#include <iostream>
#include <set>
#include <vector>
#include <deque>
#include <list>
#include "route.h"
#include "lib_record.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <algorithm>
#include <unistd.h>
#include "lp_lib.h"

using namespace std;

struct node
{
	int ID;
	int start;
	int end;
	int weight;
};

set<int> tiedot;
vector<int> dot;
int dotnum; 
int tiedotnum;
vector<int>::iterator iter;
set<int>::iterator siter;
struct node edge[5000];
int source;
int destination;
bool flag;



int dotID(int x)
{
	iter = find(dot.begin(),dot.end(),x);
	return int(iter-dot.begin());
}



//你要完成的功能总入口
void search_route(char *topo[5000], int edge_num, char *demand)
{
	
	
	source = atoi(strsep(&demand,",|"));
	destination = atoi(strsep(&demand,",|"));
	//初始化tiedot数组
	char * p = strsep(&demand,",|");
	while(NULL != p)
	{
		tiedot.insert(atoi(p));
		p = strsep(&demand,",|");
	}
	tiedotnum = tiedot.size();
	cout<<tiedotnum<<endl;
	//初始化edge数组
	for(int i=0;i<edge_num;i++)
	{
		edge[i].ID = atoi(strsep(&topo[i],","));//topo的完整性被破坏
		edge[i].start = atoi(strsep(&topo[i],","));
		edge[i].end = atoi(strsep(&topo[i],","));
		edge[i].weight = atoi(strsep(&topo[i],","));
	}
	//初始化dot数组
	for(int i=0;i < edge_num;i++)
	{
		dot.push_back(edge[i].start);
		dot.push_back(edge[i].end);
	}
	sort(dot.begin(),dot.end());
	dot.erase(unique(dot.begin(),dot.end()),dot.end());
	dotnum = dot.size();

	lprec *lp;
	int Ncol, *colno = NULL, j, ret = 0;
	REAL *row = NULL;
	Ncol = edge_num; 
	lp = make_lp(0, Ncol);
    if(lp == NULL)
		ret = 1; /* couldn't construct a new model... */
	
	set_bounds_tighter(lp, TRUE);
	
	if(ret == 0) 
	{
		
		/* create space large enough for one row */
		colno = (int *) malloc(Ncol * sizeof(*colno));
		row = (REAL *) malloc(Ncol * sizeof(*row));
		if((colno == NULL) || (row == NULL))
		  ret = 2;
  
		for(int i=1;i<=edge_num;i++)
		{
			set_binary(lp,i,true);
		}

	}
	
	if(ret == 0) //起点入度0
	{
		set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

		j = 0;
		
		for(int i=0;i<edge_num;i++)
		{
			if(edge[i].end == source)
			{
				colno[j] = i+1;
				row[j++] = 1;
			}
		}
		/* add the row to lpsolve */
		if(!add_constraintex(lp, j, row, colno, EQ, 0))
		  ret = 3;
	}
	
	if(ret == 0) //起点出度1
	{

		j = 0;
		
		for(int i=0;i<edge_num;i++)
		{
			if(edge[i].start == source&&edge[i].end != destination)
			{
				colno[j] = i+1;
				row[j++] = 1;
			}
		}
		/* add the row to lpsolve */
		if(!add_constraintex(lp, j, row, colno, EQ, 1))
		  ret = 3;
	}
	
	if(ret == 0) //终点入度1
	{

		j = 0;
		
		for(int i=0;i<edge_num;i++)
		{
			if(edge[i].end == destination&&edge[i].start != source)
			{
				colno[j] = i+1;
				row[j++] = 1;
			}
		}
		/* add the row to lpsolve */
		if(!add_constraintex(lp, j, row, colno, EQ, 1))
		  ret = 3;
	}
	
	if(ret == 0) //终点出度0
	{

		j = 0;
		
		for(int i=0;i<edge_num;i++)
		{
			if(edge[i].start == destination)
			{
				colno[j] = i+1;
				row[j++] = 1;
			}
		}
		/* add the row to lpsolve */
		if(!add_constraintex(lp, j, row, colno, EQ, 0))
		  ret = 3;
	}
	
	if(ret == 0)//必经点入度1
	{
		for(siter=tiedot.begin();siter!=tiedot.end();siter++)
		{
			j=0;
			for(int i=0;i<edge_num;i++)
			{
				if(edge[i].end == *siter)
				{
					colno[j] = i+1;
					row[j++] = 1;
				}
			}
			if(!add_constraintex(lp, j, row, colno, EQ, 1))
				ret = 3;
		}
	}
	
	if(ret == 0)//必经点出度1
	{
		for(siter=tiedot.begin();siter!=tiedot.end();siter++)
		{
			j=0;
			for(int i=0;i<edge_num;i++)
			{
				if(edge[i].start == *siter)
				{
					colno[j] = i+1;
					row[j++] = 1;
				}
			}
			if(!add_constraintex(lp, j, row, colno, EQ, 1))
				ret = 3;
		}
	}
	
	if(ret == 0) //每个点(除起始点)的出入度相等
	{

		/* construct first row (120 x + 210 y <= 15000) */
		
		
		for(int k=0;k<dotnum;k++)
		{
			if(dot[k] != source&&dot[k] != destination)
			{
				j = 0;
				for(int i=0;i<edge_num;i++)
				{
					if(edge[i].end == dot[k])
					{
						colno[j] = i+1;
						row[j++] = 1;
					}
					else if(edge[i].start == dot[k])
					{
						colno[j] = i+1;
						row[j++] = -1;
					}
				}
				
				/* add the row to lpsolve */
				if(!add_constraintex(lp, j, row, colno, EQ, 0))
				  ret = 3;
			}
			
		}
	}
	
	if(ret == 0) //每个点的入度小于等于1
	{
		for(int k=0;k<dotnum;k++)
		{
			j = 0;
			for(int i=0;i<edge_num;i++)
			{
				if(edge[i].end == dot[k])
				{
					colno[j] = i+1;
					row[j++] = 1;
				}
			}
			/* add the row to lpsolve */
			if(!add_constraintex(lp, j, row, colno, LE, 1))
			  ret = 3;
		}
	}
	
	
	if(ret == 0) //目标函数
	{
		set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

		/* set the objective function (143 x + 60 y) */
		j = 0;
		
		for(int i=0;i<edge_num;i++)
		{
			colno[j] = i+1; /* first column */
			row[j++] = edge[i].weight;
		}
		

		/* set the objective in lpsolve */
		if(!set_obj_fnex(lp, j, row, colno))
		  ret = 4;
	}
	
	vector<int> tempresult;
	set_minim(lp);
	if(dotnum>=550)
	{
		set_break_at_first(lp,TRUE);
	}
	set_verbose(lp, IMPORTANT);
	do
	{
		flag = false;
		
		if(ret == 0||ret == 5) 
		{
			ret = solve(lp);
			if(ret == OPTIMAL)
			  ret = 0;
			else
			  ret = 5;
		}
		
		
		if(ret == 0||ret == 5) 
		{
			printf("Objective value: %f\n", get_objective(lp));

			set<int> s;
			set<int> stemp;
			set<int> circledot;
			set<int>::iterator it;
			
			get_variables(lp, row);
			for(j = 0; j < Ncol; j++)
			{
				if(row[j] == 1)
				{
					s.insert(j);
				}
			}
			
			int tempid = source;
			while(!s.empty())
			{
				for(it=s.begin();it!=s.end();it++)
				{
					if(edge[*it].start == tempid)
					{
						tempresult.push_back(edge[*it].ID);
						tempid = edge[*it].end;
						s.erase(*it);
						break;
					}
					
					if(!s.empty()&&it == (--s.end()))
					{
						for(it=s.begin();it!=s.end();it++)
						{
							circledot.insert(edge[*it].start);
						}
						
						circledot.insert(edge[*(--it)].end);
						stemp = s;
						s.clear();
						tempresult.clear();
						
						flag = true;
						break;
					}
				}
			}
			//取子环
			if(!circledot.empty())
			{
				
				do
				{
					
					set<int> intersection;
					set<int> diff;
					set_intersection(circledot.begin(),circledot.end(),tiedot.begin(),tiedot.end(),
									inserter(intersection,intersection.begin()));
					set<int> subcircledot;
					tempid = *intersection.begin();
					
					do
					{
						for(it=stemp.begin();it!=stemp.end();it++)
						{
							if(edge[*it].start == tempid)
							{
								subcircledot.insert(tempid);
								tempid = edge[*it].end;
								stemp.erase(*it);
								break;
							}
						}
					}while(tempid != *intersection.begin());
					
					
					
					j=0;
					for(int i=0;i<edge_num;i++)
					{
						if(subcircledot.find(edge[i].start)!=subcircledot.end()&&
						   subcircledot.find(edge[i].end)!=subcircledot.end())
						{
							colno[j] = i+1;
							row[j++] = 1;
						}
					}
					
					if(!add_constraintex(lp, j, row, colno, LE, subcircledot.size()-1))
					  ret = 3;

					set_difference(circledot.begin(),circledot.end(),subcircledot.begin(),subcircledot.end(),
									inserter(diff,diff.begin()));
					
					circledot = diff;
					
					
				}while(!circledot.empty());
				
				
			}
		}

	}while(flag);
	
	for(int i=0;i<tempresult.size();i++)
	{
		record_result(tempresult[i]);
	}
	
	/* free allocated memory */
	if(row != NULL)
	free(row);
	if(colno != NULL)
	free(colno);

	if(lp != NULL) 
	{
	/* clean up such that all used memory by lpsolve is freed */
	delete_lp(lp);
	}
 
}



