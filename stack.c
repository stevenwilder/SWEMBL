#include "SWEMBL.h"

void initialize(stack *stk)
{
  stk -> cnt = 0;
  stk -> head = NULL;
  stk -> tail = NULL;
}

void unshift(long int start, long int end, double count, stack *stk)
{
  elem *p;
  p = malloc(sizeof(elem));
  p -> d[0] = start;
  p -> d[1] = end;
  p -> s = count;
  p -> next = stk -> head;
  stk -> head = p;
  if(stk -> cnt == 0) { stk -> tail = p; }
  stk -> cnt++;
}

int cnt(stack *stk)
{
  return(stk -> cnt);
}



void push(long int start, long int end, double count, stack *stk)
{
  /*if(stk -> tail && stk -> tail -> d[0] == start && stk -> tail -> d[1] == end)
    {
      stk -> tail -> s += count;
    }

  else
  {*/
      elem *p;
      p = malloc(sizeof(elem));
      p -> d[0] = start;
      p -> d[1] = end;
      p -> s = count;
      //if(stk -> tail)
      if(stk -> cnt > 0)
	{ stk -> tail -> next = p; }
      else
	{ stk -> head = p;}
      stk -> tail = p;
      stk -> cnt++;
      // }
}

struct startendcount shift(stack *stk)
{
  elem *p;
  struct startendcount d;
  p = stk -> head;
  d.startend[0] = p -> d[0];
  d.startend[1] = p -> d[1];
  d.count = p -> s;
  //stk -> head = stk -> head -> next;
  stk -> cnt--;
  //printf("shc:%d\n", stk->cnt);
  //elem *q;
  //q = NULL;
  if(!(stk -> cnt)) {stk -> head -> next = NULL; stk -> tail = NULL;}
  //if(!(stk -> cnt)) { stk -> tail = NULL ; }
  else {stk -> head = stk -> head -> next;}
  free(p);
  //if(stk -> cnt) { free(p); }
  return d;
}

void add_to_position(double n, int posn, stack *stk)
{
  elem *p;
  switch(posn)
    {
    case '0': p = stk->head; break;
    default: 
      {
	if(posn >= cnt(stk)) { return; }
	else{ if( posn == (cnt(stk) - 1)) {p =stk->tail;}
	else
	  {
	    p = stk->head;
	    int i;
	    for(i =0; i < posn; i++)
	      { 
		//printf("add\t%d\t%d\t%d\n", i, posn, cnt(stk));
		p = p->next; 
	      }
	  }
	}
      }
    }
  
  /*for(i =0; i < posn; i++)
     { 
      // printf("add\t%d\t%d\t%d\n", i, posn, cnt(stk));
       p = p->next; 
     }
  */
  //(p -> d[start_or_end]) += n;
  (p -> s) += n; 
}

void add_to_or_insert_start_position(long int start, long int end, double n, stack *stk)
{
  if(empty(stk))
    {
      //printf("Emp\t%ld\t%ld\t%f\n", start, end, n);
      push(start, end, n, stk);
      return;
    }

  elem *p;
  p = stk -> tail;

  //printf("AISPE %ld\t%ld\t%ld\t%ld\t%f\n", start, end, stk->head->d[0], p->d[0], n); 
  if(start != p->d[0])
    {
      if(start > p->d[0] )
	{
	  //printf("Tail\t%ld\t%ld\t%f\n", start, end, n);
	  push(start, end, n, stk);
	  return;
	}
      else
	{
	  p = stk -> head;
	  if(start < p->d[0])
	    {	      
	      //printf("Hd\t%ld\t%ld\t%f\t%ld\n", start, end, n, stk->head->d[0]);
	      unshift(start, end, n, stk);
	      return;
	    }	      
	  int i = 1;
	  while(i < (stk->cnt) && p->next->d[0] <= start)
	    {
	      p = p->next;
	      i++;
	    }
	  if(p->d[0] < start)
	    {
	      elem *q;
	      q = malloc(sizeof(elem));
	      q -> s = n;
	      q -> d[0] = start;
	      q -> d[1] = end;
	      q -> next = p -> next;
	      p -> next = q;
	      stk -> cnt++;
	      //printf("Ins\t%ld\t%ld\t%f\n", start, end, n);
	      return;
	    }	      
	}
    }
  (p -> s) += n;
  return;
}

  

long int return_position(int posn, stack *stk)
{
  elem *p;
  switch(posn)
    {
    case '0': p = stk->head; return(p->d[0]); break;
    default: 
      {
	if(posn >= cnt(stk)) { return 0; }
	else{ if( posn == (cnt(stk) - 1)) {p =stk->tail; return p->d[0];}
	else
	  {
	    p = stk->head;
	    int i;
	    for(i =0; i < posn; i++)
	      { 
		//printf("ret\t%d\t%d\t%d\n", i, posn, cnt(stk));
		p = p->next; 
	      }
	    return p->d[0];
	  }
	}
      }
    }
}

struct startendcount head(stack *stk)
{
  struct startendcount d;
  d.startend[0] = stk -> head -> d[0];
  d.startend[1] = stk -> head -> d[1];
  d.count = stk -> head -> s;
  return (d);
}

struct startendcount tail(stack *stk)
{
  struct startendcount d;
  d.startend[0] = stk -> tail -> d[0];
  d.startend[1] = stk -> tail -> d[1];
  d.count = stk -> tail -> s;
  return (d);
}

boolean empty(const stack *stk)
{
  return ((boolean) (stk -> cnt == 0));
}

void assign(stack *stk, stack *istk)
{
  *stk = *istk;
}

struct countend cntequalhead(long int pos, stack *stk)
{
  int n = 0;
  double count = 0;
  long int end = 0;
  if(stk->cnt != 0)
    {
      elem *q;
      q = stk -> head;
      if(q->d[0] == pos)
	{
	  n = 1;
	  count = q->s;
	  end = q->d[1];
	  while(n < stk->cnt && (q -> next -> d[0] == pos))
	    {
	      n++;
	      count += q->s;
	      if(q -> next -> d[1] > end){ end = q -> next -> d[1] ; }
	      if(q->next)
		{
		  q = q->next;		  
		}
	      else{break;}
	    }
	}
    }
  struct countend ce;
  ce.count = count;
  ce.end = end;
  return(ce);
}

// *** Sort stk by highest end first (head) then by lowest start ***
void sortunshift(long int start, long int end, double count, stack *stk)
{
  elem *p;
  p = malloc(sizeof(elem));
  p -> d[0] = start;
  p -> d[1] = end;
  p -> s = count;

  //if(stk->cnt){printf("sus: %ld\t%ld\t%ld\t%ld\t%f\n", start, end, stk->head->d[0], stk->tail->d[0], count);}

  // ***Read goes to head of savedpos ***
  if((stk->cnt==0) || (end > stk->head->d[1]) || (end == stk->head->d[1] && start <= stk->head->d[0]))
    {
      p -> next = stk -> head;
      stk -> head = p;
      if(stk->cnt == 0) { stk->tail = p; p->next=NULL;}
      //printf("a:%ld\t%ld\n", start, end);
    }
 else
    {
      // ***Read goes to tail of savedpos***
      //if(stk->tail && (end < stk->tail->d[1] || (end == stk->tail->d[1] && start <= stk->tail->d[0])))
      //if(stk->tail && (end < stk->tail->d[1] || (end == stk->tail->d[1] && start >= stk->tail->d[0])))
      if(stk->tail && (end < stk->tail->d[1] || (end == stk->tail->d[1] && start >= stk->tail->d[0])))
	{
	  if(stk -> tail)
	    { stk -> tail -> next = p; }
	  else
	    { stk -> head = p;}
	  stk -> tail = p;
	  //printf("t:%ld\t%ld\n", start, end);
	}

      //***Read inserted into non-terminal position of savedpos***
      else
	{
	  elem *q;
	  q = stk->head;
	  while(q->next && end <= q->next->d[1] && start > q->next->d[0])
	    {
	      //printf("r:%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n", q->next->d[0],q->next->d[1], start, end,stk->head->d[0],stk->head->d[1],stk->tail->d[0],stk->tail->d[1]);
	    ////printf("r:%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n", q->next->d[0],q->next->d[1], start, end,stk->head->d[0],stk->head->d[1]);
	      q = q->next;
	    }
	  p -> next = q -> next;
	  q -> next = p; 
	}
    }
  stk->cnt++;
}

double median(stack *stk)
{
  if(empty(stk))
    { printf("Median of empty stack taken!\n"); exit(1); }

  double midlength = ((double)cnt(stk) / 2);
  //printf("C %d\tM %f\n", cnt(stk), midlength);
  double med;
  double i;
  
  for(i = 1; i < midlength; i++)
    {
      shift(stk);
    }

  if((double)floor(midlength) == midlength) // ***EVEN***
    {
      //med = shift(stk).count/2;
      //med += shift(stk).count/2;  
      med = (double) shift(stk).startend[0]/2;
      med += (double) shift(stk).startend[0]/2;
    }
  else
    {
      //med = shift(stk).count;
      med = (double) shift(stk).startend[0];
    }
  return med;
}

void initializeregion(regionstack *stk)
{
  stk -> cnt = 0;
  stk -> head = NULL;
  stk -> tail = NULL;
}

int regioncnt(regionstack *stk)
{
  return(stk -> cnt);
}

void unshiftregion(struct region r, regionstack *stk)
{
  regionelem *p;

  p = malloc(sizeof(regionelem));
  p -> r = r;
  p -> next = stk -> head;
  stk -> head = p;
  stk -> cnt++;
}


struct region shiftregion(regionstack *stk)
{
  struct region r;
  regionelem *p;
  
  r = stk -> head -> r;
  p = stk -> head;
  stk -> head = stk -> head -> next;
  stk -> cnt--;
  if(!(stk -> cnt)) { stk -> tail = NULL;}
  free(p);
  return r;
}

struct region headregion(regionstack *stk)
{
  return (stk -> head -> r);
}

boolean emptyregion(const regionstack *stk)
{
  return ((boolean) (stk -> cnt == 0));
}


void initializechr(chrstack *stk)
{
  stk -> cnt = 0;
  stk -> head = NULL;
  stk -> tail = NULL;
}

int chrcnt(chrstack *stk)
{
  return(stk -> cnt);
}

void pushchr(long int start, long int end, double count, chrstack *cstk)
{
  stack *posstk = lastchr(cstk);
  elem *p;

  p = malloc(sizeof(elem));
  p -> d[0] = start;
  p -> d[1] = end;
  p -> s = count;

  if(posstk -> tail)
    { posstk -> tail -> next = p; }
  else
    { posstk -> head = p;}
  posstk -> tail = p;
  posstk -> cnt++;
}

void addchr(char chr[50], chrstack *cstk)
{
  chrelem *p;

  p = malloc(sizeof(chrelem));
  strcpy(p -> chr, chr);
  //printf("AC\n");

  stack *s; 
  s = malloc(sizeof(elem));
  s -> cnt = 0;  
  //printf("AC1\n");
  s -> head = NULL;
  s -> tail = NULL;
  
  //printf("AC2\n");

  p->stk = s;
  //printf("AC3\n");
  if(cstk -> tail)
    { cstk -> tail -> next = p; }
  else
    { cstk -> head = p;}
  cstk -> tail = p;
  cstk -> cnt++;
  //printf("AC5\n");
}

struct startendcount shiftchr(chrstack *cstk)
{
  stack *posstk = firstchr(cstk);
  //printf("SC1\n");
  elem *p;
  struct startendcount d;
  p = posstk -> head;
  d.startend[0] = p -> d[0];
  d.startend[1] = p -> d[1];
  d.count = p -> s;
  posstk -> head = posstk -> head -> next;
  //printf("SC2\t%ld\n",d.startend[0]);
  posstk -> cnt--;
  if(!(posstk -> cnt)) { posstk -> tail = NULL;}
  
  //printf("SC3\n");
  strcpy(prev.compchr, firstchrelem(cstk)->chr);
  prev.compstart = d.startend[0];
  prev.compend = d.startend[1];
  free(p);
  
  
  //if(posstk -> cnt == 0)
  // { removechr(cstk); }
  //printf("SC4\n");

  return d;
}

struct stack *removechr(chrstack *chrstk)
{
  struct stack *posstk;
  chrelem *p;
  
  posstk = chrstk -> head -> stk;
  p = chrstk -> head;
  chrstk -> head = chrstk -> head -> next;
  chrstk -> cnt--;
  if(!(chrstk -> cnt)) { chrstk -> tail = NULL;}
  free(p);
  return posstk;
}

struct startendcount headchr(chrstack *cstk)
{
  stack *posstk = firstchr(cstk);
  struct startendcount d;
  d.startend[0] = posstk -> head -> d[0];
  d.startend[1] = posstk -> head -> d[1]; 
  d.count = posstk -> head -> s;
  return(d);
}

struct startendcount nextheadchr(chrstack *cstk)
{
  stack *posstk = firstchr(cstk);
  struct startendcount d;
  d.startend[0] = posstk -> head -> next -> d[0];
  d.startend[1] = posstk -> head -> next -> d[1];
  d.count = posstk -> head -> next -> s;
  return(d);
}


struct stack *firstchr(chrstack *cstk)
{ 
  if(!emptychr(cstk))
    {return (cstk -> head -> stk);}
  else
    {
      struct stack *p;
      p=malloc(sizeof(elem));
      p -> cnt = 0;
      p -> head = NULL;
      p -> tail = NULL;
      return(p);
    }
}


struct chrelem *firstchrelem(chrstack *cstk)
{
  return (cstk -> head);
}

struct chrelem *secondchrelem(chrstack *cstk)
{
  if(chrcnt(cstk) > 1)
    {
      return(cstk -> head -> next); 
    }
  else
    {
      struct chrelem *c;
      c=malloc(sizeof(chrelem));
      strcpy(c->chr,"");
      return(c);
    }
}

struct stack *lastchr(chrstack *cstk)
{ 
  if(!emptychr(cstk))
    {return (cstk -> tail -> stk);}
  else
    {
      struct stack *p;
      p=malloc(sizeof(elem));
      p -> cnt = 0;
      p -> head = NULL;
      p -> tail = NULL;
      return(p);
    }
}

struct chrelem *lastchrelem(chrstack *cstk)
{
  if(!emptychr(cstk))
    {return (cstk -> tail);}
  else
    {
      struct chrelem *c;
      c=malloc(sizeof(chrelem));
      strcpy(c->chr,"");
      return(c);
    }
}


boolean emptychr(const chrstack *stk)
{
  return ((boolean) (stk -> cnt == 0));
}


int seenchr(char chr[], struct chrstack *cstk)
{
  int i;
  chrelem *c;
  c = cstk -> head;

  for(i=0; i< chrcnt(cstk); i++)
    {
      if(!strcmp(chr, c->chr))
	{return(1);}
      c = c -> next;
    }
  return(0);
}
