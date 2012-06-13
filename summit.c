#include "SWEMBL.h"

void push_peakheight(int dir, long int start, long int end, double score, struct region *curr, struct region *possend)
{
  //int i;
  //printf("pph %d %ld %ld %f\n", dir, start, end, score);

  if(!empty(&peakheights))
    {     
      calculate_peak_summit(curr, head(&peakheights).startend[0], dir * start - 1, possend);
    }

  if(!empty(&peakheights))
    {
      //add_to_position(score, 0, &peakheights);
      add_to_or_insert_start_position(dir * start, dir * end, score, &peakheights);
      //printf("aisp %f\t%ld\t%ld\n", score, dir * start, head(&peakheights).startend[0]);
    }
  else
    {
      push(dir * start, dir * end, score, &peakheights);
      //printf("ps. %f\t%ld\n", head(&peakheights).count, head(&peakheights).startend[0]);
    }

 
  /*
    for(i = dir * start + cnt(&peakheights); i <= dir * end; i++)
     {
       // printf("pushi %d\n", i);
       push(dir*i, 0, 0, &peakheights);
     }
  */
  
  //add_to_position(-1 * score, (end - start) * dir , &peakheights);
  
  add_to_or_insert_start_position(dir * end + 1, dir * end + 1, -1 * score, &peakheights);
  //printf("aispe %f\t%ld\t%ld\n", -1 *score, dir * end + 1, head(&peakheights).startend[0]);
} 

	     
void calculate_peak_summit(struct region *curr, long int begin, long int finish, struct region *possend)
{ 
  int i;
  //begin *= dir;
  //finish *= dir;
  //dir = 1;
  //for(i = dir * head(&peakheights).startend[0]; i <= dir * tail(&peakheights).startend[0]; i++)

  //for(i = ((dir *head(&peakheights).startend[0]) < (dir * begin)) ? (dir *head(&peakheights).startend[0]) : (dir * begin)    ; i <= (dir * finish); i++)
  for(i = (head(&peakheights).startend[0] < begin) ? head(&peakheights).startend[0] : begin ; i <= finish; i++)
    {    
      if(cnt(&peakheights) == 0)
	{break;}

            

      
      if(i == head(&peakheights).startend[0])    //NEW
	{ 
	  (*curr).peakheight += shift(&peakheights).count;      
	  (*possend).peakheight = (*curr).peakheight;
	  //printf("%d\t%f\n", i, (*curr).peakheight);
	}

      // printf("i %d\tbegin %ld\tfinish %ld\ts %f\tp %f\tm %f\tpossend %f\t%ld\tb %ld\n", i, ((dir *head(&peakheights).startend[0]) < (dir * begin)) ? (dir *head(&peakheights).startend[0]) : (dir * begin), (dir * finish), head(&peakheights).count, (*curr).peakheight, (*curr).peakheightmax, (*possend).peakheightmax, head(&peakheights).startend[0], dir*begin);

      //printf("i %d\tbegin %ld\tfinish %ld\ts %f\tp %f\tm %f\tpossend %f\t%ld\tb %ld\n", i, (head(&peakheights).startend[0] < begin) ? head(&peakheights).startend[0] : begin, finish, head(&peakheights).count, (*curr).peakheight, (*curr).peakheightmax, (*possend).peakheightmax, head(&peakheights).startend[0], begin);


      if((*curr).peakheight > (*curr).peakheightmax || empty(&peaksummits))
	{
	  (*curr).peakheightmax = (*curr).peakheight;
	  (*possend).peakheightmax = (*curr).peakheight;
	  while(!empty(&peaksummits)) { shift(&peaksummits); }
	  push(i, i, (*curr).peakheightmax, &peaksummits);
	}

      else
	{ 
	  if((*curr).peakheight == (*curr).peakheightmax)
	    {   
	      push(i, i, (*curr).peakheightmax, &peaksummits);
	    }
	}
    }
}
