
#include "gft_graph.h"

namespace gft{
  namespace Graph{


    Graph *Create(int nnodes, int outdegree, float *Wnodes){
      Graph *graph = NULL;
      int i;
      graph = (Graph *) calloc(1, sizeof(Graph));
      if(graph == NULL){
	gft::Error((char *)MSG1,(char *)"Graph::Create");
      }
      graph->nnodes = nnodes;
      graph->Wnodes = Wnodes;
      graph->nodes = (GraphNode *) calloc(nnodes, sizeof(GraphNode));
      if(graph->nodes == NULL){
	gft::Error((char *)MSG1,(char *)"Graph::Create");
      }
      for(i = 0; i < nnodes; i++){
	graph->nodes[i].outdegree = 0;
	graph->nodes[i].arraysize = outdegree;
	graph->nodes[i].adjList = gft::AllocIntArray(outdegree);
	graph->nodes[i].Warcs = gft::AllocFloatArray(outdegree);
      }
      return graph;
    }

    
    Graph *Clone(Graph *graph){
      Graph *clone;
      int p,q,i;
      float w;
      clone = Create(graph->nnodes, 10, NULL);
      if(graph->Wnodes != NULL){
	clone->Wnodes = (float *)malloc(sizeof(float)*graph->nnodes);
	for(p = 0; p < graph->nnodes; p++){
	  clone->Wnodes[p] = graph->Wnodes[p];
	}
      }
      for(p = 0; p < graph->nnodes; p++){
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  w = graph->nodes[p].Warcs[i];
	  AddDirectedEdge(clone, p, q, w);
	}
      }
      return clone;
    }


    Graph *Transpose(Graph *graph){
      Graph *transp;
      int p,q,i;
      float w;
      transp = Create(graph->nnodes, 10, NULL);
      if(graph->Wnodes != NULL){
	transp->Wnodes = (float *)malloc(sizeof(float)*graph->nnodes);
	for(p = 0; p < graph->nnodes; p++){
	  transp->Wnodes[p] = graph->Wnodes[p];
	}
      }
      for(p = 0; p < graph->nnodes; p++){
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  w = graph->nodes[p].Warcs[i];
	  AddDirectedEdge(transp, q, p, w);
	}
      }
      return transp;
    }    
    

    void   Destroy(Graph **graph){
      Graph *aux;
      int i;
      if(graph != NULL){
	aux = *graph;
	if(aux != NULL){
	  if(aux->nodes != NULL){
	    for(i = 0; i < aux->nnodes; i++){
	      if(aux->nodes[i].adjList != NULL)
		free(aux->nodes[i].adjList);
	      if(aux->nodes[i].Warcs != NULL)
		free(aux->nodes[i].Warcs);
	    }
	    free(aux->nodes);
	  }
	  free(aux);
	  *graph = NULL;
	}
      }
    }



    float GetArcWeight(Graph *graph, int src, int dest){
      GraphNode *s;
      s = &graph->nodes[src];
      int n,i;
      n = s->outdegree;
      for(i = 0; i < n; i++){
	if(s->adjList[i] == dest){
		//printf("s->Warcs %f\n",s->Warcs[i]);
	  	return s->Warcs[i];
	}
      }
      gft::Warning((char *)"Invalid node", (char *)"Graph::GetArcWeight");
      return -1.0;
    }
    
    
    void AddEdge(Graph *graph, int src, int dest, float w){
      AddDirectedEdge(graph, src,  dest, w);
      AddDirectedEdge(graph, dest, src,  w);  
    }
    
    
    void AddDirectedEdge(Graph *graph, int src, int dest, float w){
      GraphNode *s;
      s = &graph->nodes[src];
      s->outdegree++;
      if(s->outdegree > s->arraysize){
	s->adjList = (int *)realloc(s->adjList, sizeof(int)*s->outdegree);
	if(s->adjList == NULL)
	  gft::Error((char *)MSG1,(char *)"Graph::AddDirectedEdge");
	s->Warcs = (float *)realloc(s->Warcs, sizeof(float)*s->outdegree);
	if(s->Warcs == NULL)
	  gft::Error((char *)MSG1,(char *)"Graph::AddDirectedEdge");
	s->arraysize = s->outdegree;
      }
      s->adjList[s->outdegree-1] = dest;
      s->Warcs[s->outdegree-1] = w;
      //printf("W_graph %f\n",s->Warcs[s->outdegree-1]);
    }
    
    
    void UpdateEdge(Graph *graph, int src, int dest, float w){
      UpdateDirectedEdge(graph, src,  dest, w);
      UpdateDirectedEdge(graph, dest, src,  w);  
    }
    

    void UpdateDirectedEdge(Graph *graph, int src, int dest, float w){
      GraphNode *s;
      int i, find = false;
      s = &graph->nodes[src];

      for(i = 0; i < s->outdegree; i++){
	if(dest == s->adjList[i]){
	  s->Warcs[i] = w;
	  find = true;
	}
      }
      if(!find)
	AddDirectedEdge(graph, src, dest, w);
    }


    void UpdateDirectedEdgeIfHigher(Graph *graph, int src, int dest, float w){
      GraphNode *s;
      int i, find = false;
      s = &graph->nodes[src];

      for(i = 0; i < s->outdegree; i++){
	if(dest == s->adjList[i]){
	  if(s->Warcs[i] < w)
	    s->Warcs[i] = w;
	  find = true;
	  break;
	}
      }
      if(!find)
	AddDirectedEdge(graph, src, dest, w);
    }


    float GetMaximumArc(Graph *graph){
      float w,Wmax = FLT_MIN;
      int p,i;
      for(p = 0; p < graph->nnodes; p++){
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  w = graph->nodes[p].Warcs[i];
	  if(w > Wmax)
	    Wmax = w;
	}
      }
      return Wmax;
    }
    

    int Tarjan_aux(Graph *graph, int *index, int *lowlink, int *label,
		   int i, int id, int *lb){
      //pilhas usadas para simular a recursao:
      gft::Stack::Stack *S_i = NULL;
      gft::Stack::Stack *S_j = NULL;
      //pilha do algoritmo de Tarjan:
      gft::Stack::Stack *S = NULL;
      int *IsInS = NULL;
      int j, t, k, tmp, flag = 0;
      enum action {f_in, f_out, f_end} act;

      IsInS = gft::AllocIntArray(graph->nnodes);
      S_i = gft::Stack::Create(graph->nnodes);
      S_j = gft::Stack::Create(graph->nnodes);
      S   = gft::Stack::Create(graph->nnodes);
      act = f_in;
      //---------------------------
      j = 0;
      gft::Stack::Push(S_i, i);
      gft::Stack::Push(S_j, j);
      index[i] = id;
      lowlink[i] = id;
      id++;
      gft::Stack::Push(S, i);
      IsInS[i] = 1;

      //--------------------------
      do{
	switch(act){
	case f_in:
	  if(j < graph->nodes[i].outdegree){
	    k = graph->nodes[i].adjList[j];
	    j = 0;
	    i = k;
	  }
	  else{
	    act = f_out;
	    break;
	  }

	  if(index[i] != NIL){
	    flag = 1;
	    if(IsInS[i]){
	      //node.lowlink = min(node.lowlink, n.index);
	      tmp = gft::Stack::Pop(S_i);
	      if(lowlink[tmp] > index[i])
		lowlink[tmp] = index[i];
	      gft::Stack::Push(S_i, tmp);
	    }
	    act = f_out;
	  }
	  else{
	    index[i] = id;
	    lowlink[i] = id;
	    id++;
	    gft::Stack::Push(S, i);
	    IsInS[i] = 1;
	    gft::Stack::Push(S_i, i);
	    gft::Stack::Push(S_j, j);
	  }
	  break;
	case f_out:
	  if(gft::Stack::IsEmpty(S_i))
	    act = f_end;
	  else{
	    tmp = lowlink[i];
	    i = gft::Stack::Pop(S_i);
	    j = gft::Stack::Pop(S_j);

	    if(lowlink[i] > tmp && !flag)
	      lowlink[i] = tmp;

	    flag = 0;
	    j++;
	    if(j < graph->nodes[i].outdegree){
	      gft::Stack::Push(S_i, i);
	      gft::Stack::Push(S_j, j);
	      act = f_in;
	    }
	    else{
	      //if we are in the root of the component
	      if(lowlink[i] == index[i]){
		do{
		  if(gft::Stack::IsEmpty(S))
		    break;
		  t = gft::Stack::Pop(S);
		  IsInS[t] = 0;
		  label[t] = *lb;
		}while(t != i);
		(*lb)++;
	      }
	    }
	  }
	  break;
	case f_end:
	  break;
	}
      }while(act != f_end);

      gft::FreeIntArray(&IsInS);
      gft::Stack::Destroy(&S_i);
      gft::Stack::Destroy(&S_j);
      gft::Stack::Destroy(&S);
      return id;
    }

    
    /*Tarjan's algorithm is a procedure for finding strongly connected components of a directed graph.*/
    int *Tarjan(Graph *graph){
      int *index = NULL;
      int *lowlink = NULL;
      int *label = NULL;
      int i,id = 0,lb = 1;

      index   = gft::AllocIntArray(graph->nnodes);
      lowlink = gft::AllocIntArray(graph->nnodes);
      label   = gft::AllocIntArray(graph->nnodes);
      for(i = 0; i < graph->nnodes; i++)
	index[i] = NIL;
      	
      for(i = 0; i < graph->nnodes; i++){
	if(index[i] != NIL) continue;
	id = Tarjan_aux(graph, index, lowlink, label, i, id, &lb);
      }

      gft::FreeIntArray(&index);
      gft::FreeIntArray(&lowlink);
      return label;
    }


    int *Tarjan(Graph *graph, int p){
      int *index = NULL;
      int *lowlink = NULL;
      int *label = NULL;
      int id = 0,lb = 1, i;

      index   = gft::AllocIntArray(graph->nnodes);
      lowlink = gft::AllocIntArray(graph->nnodes);
      label   = gft::AllocIntArray(graph->nnodes);
      for(i = 0; i < graph->nnodes; i++)
	index[i] = NIL;

      id = Tarjan_aux(graph, index, lowlink, label, p, id, &lb);

      lb = label[p];
      for(i = 0; i < graph->nnodes; i++){
	if(label[i] != lb)
	  label[i] = 0;
	else
	  label[i] = 1;
      }
      
      gft::FreeIntArray(&index);
      gft::FreeIntArray(&lowlink);
      return label;
    }


    int *Tarjan(Graph *graph, int *V, int n){
      int *index = NULL;
      int *lowlink = NULL;
      int *label = NULL;
      int i,id = 0,lb = 1;

      index   = gft::AllocIntArray(graph->nnodes);
      lowlink = gft::AllocIntArray(graph->nnodes);
      label   = gft::AllocIntArray(graph->nnodes);
      for(i = 0; i < graph->nnodes; i++)
	index[i] = NIL;
      	
      for(i = 0; i < n; i++){
	if(index[V[i]] != NIL) continue;
	id = Tarjan_aux(graph, index, lowlink, label, V[i], id, &lb);
      }

      //----------------------------
      int *lb_map;
      int Lmax = 0,l;
      for(i = 0; i < graph->nnodes; i++){
	if(label[i] == NIL)
	  label[i] = 0;
	if(label[i] > Lmax)
	  Lmax = label[i];
      }
      lb_map = gft::AllocIntArray(Lmax+1);

      for(i = 0; i < n; i++)
	lb_map[label[V[i]]] = 1;
      
      l = 0;
      for(i = 1; i < Lmax+1; i++){
	if(lb_map[i] > 0){
	  l++;
	  lb_map[i] = l;
	}
      }

      for(i = 0; i < graph->nnodes; i++)
	label[i] = lb_map[label[i]];

      gft::FreeIntArray(&lb_map);
      //----------------------------      

      gft::FreeIntArray(&index);
      gft::FreeIntArray(&lowlink);
      return label;
    }

    
  } /*end Graph namespace*/
} /*end gft namespace*/

