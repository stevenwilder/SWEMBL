#include "SWEMBL.h"

long int set_negpos(struct readinfo *read, int *endfile, int seqlength, int fraglength)
  {
    int changechr = 0;
    long int negpos = (*read).pos + seqlength - fraglength;
    if(*endfile && !empty(&pvepos))
      { 
	negpos = tail(&pvepos).startend[0];
	changechr = seqlength;
	//printf("End of file reached.\n");
      }
    if((*read).changechr)
      {	 
	if(!empty(&pvepos))
	  { negpos = tail(&pvepos).startend[0];}
	changechr = seqlength;
	//printf("ChC:%ld\n", negpos);
	//read.changechr = 0;
      }
    (*read).changechr = changechr;
    return(negpos);
  }

void set_fragpos(struct readinfo read, long int *fragpos, long int *fragend, long int negpos, long int fraglength)
   {
     /*
     if(!empty(&pvepos) && head(&pvepos).startend[0] < negpos)
	{ 
	  *fragpos = head(&pvepos).startend[0]; 
	  *fragend = head(&pvepos).startend[1];
	}
      else
	{ 
	  *fragpos = negpos; 
	  *fragend = negpos + fraglength - 1;
	}
     */

     long int fpos, fend;
     //if(!empty(&pvepos) && head(&pvepos).startend[0] < negpos)
     if(!empty(&pvepos) && head(&pvepos).startend[0] <= negpos)
	{ 
	  fpos = head(&pvepos).startend[0]; 
	  fend = head(&pvepos).startend[1];
	}
      else
	{ 
	  fpos = negpos; 
	  fend = negpos + fraglength - 1;
	}

    
     if(read.changechr && *fragpos <= (read.pos + read.changechr - fraglength) && fpos > (read.pos + read.changechr - fraglength) && read.negcount > 0)
       {
	 fpos = read.pos + read.changechr - fraglength;
	 fend = read.pos + fraglength - 1;
       }
     *fragpos = fpos;
     *fragend = fend;
   }


void end_chr(struct readinfo read, struct region *curr, struct param par, int fragpos, int fragend, struct region *lastregion)
   {
     if(!par.quietmode){printf("%s\n", read.chr);}
     if(strcmp(read.prevchr,""))
       {
	 //genomelength++;
	 //pvestrandregion.hightot=0;
	 //strcpy(read.chr, read.prevchr);
	 store_results(read, *curr, par, lastregion); 
       }
     
     //else { strcpy(read.prevchr, read.chr) ; }
     *curr = new_region(read, par, fragpos, fragend, 1);
   }



void update_total(int dir, struct region *curr, struct region *possend, long int negpos, long int fragpos, long int fragend, long int prevfragpos, struct param par, struct readinfo read, struct region *lastregion)
   {
     int gap = dir * (fragpos - prevfragpos);
		  
     //printf("A4:%d\t%ld\n",gap,prevfragpos);
     //genomelength = genomelength + gap;

     double newtot = (*curr).total;

     if(read.count > 0)
       {
	 newtot = (gap < par.pen_inc) ?
	   (*curr).total - par.bg - par.posbg * gap  :
	   (*curr).total - par.bg - par.posbg * par.pen_inc - par.longbg * (gap - par.pen_inc) ;

	 
	 (*curr).print_total = (gap < par.pen_inc) ?
	   (*curr).print_total - par.bg - par.posbg * gap  :
	   (*curr).print_total - par.bg - par.posbg * par.pen_inc - par.longbg * (gap - par.pen_inc) ;
       }
     
     //printf("A5\tnewtot:%f\n",newtot);
	      
     //printf("A6\n");
     if(newtot <= 0.00001)  //End of region
       {
	 //printf("A6a\n");
	 if(dir == 1) 
	   { 
	     store_results(read, *possend, par, lastregion); 
	     *curr = new_region(read, par, fragpos, fragend, dir);
	     *possend = *curr;
	   }
	 else{ (*curr).total = newtot; return; }
       }
     
     else 
       {
	 newtot = newtot + read.count - read.refcount * par.refpen;
	 (*curr).print_total += read.count - read.refcount * par.refpen;
	 //printf("A6\tnewtot:%f\t%ld\t%ld\t%ld\n",newtot,fragpos,fragend,negpos);
	 (*curr).count += read.count;
	 (*curr).refcount += read.refcount;
	 if(read.count > 0) 
	   { 
	     (*curr).pos++; 
	     //(*curr).end = fragend;
	     if(dir * fragend > dir * (*curr).end) { (*curr).end = fragend; }
	     if((negpos != fragpos || (read.count > read.negcount)) && dir * fragpos < dir * (*curr).seen_start) 
	       { 		 
		 (*curr).seen_start = fragpos; 
		 //printf("A\n"); 
	       }
	     //if(!par.paired && dir * fragpos + read.pairlength - 1 > dir * (*curr).seen_end) 
	     // { 
	     //	 (*curr).seen_end = fragpos + dir*(read.pairlength - 1); 
	     // //printf("B\n");
	      // } 
	     if(!par.paired && dir * fragpos + par.seqlength - 1 > dir * (*curr).seen_end) 
	      { 
		(*curr).seen_end = fragpos + dir*((long int)par.seqlength - 1); 
		//printf("B\t%ld\n", (*curr).seen_end);
	      } 	     
	     if(par.paired && dir * fragend > dir * (*curr).seen_end) 
	       { 
		 (*curr).seen_end = fragend ; 
		 //printf("Ba\t%ld\n", (*curr).seen_end); 
	       }
	     //if(negpos == fragpos && dir * fragend > dir * (*curr).seen_end) { (*curr).seen_end = fragend; printf("C\n");}
	     //if((dir == -1 || empty(&pvepos) || head(&pvepos).startend[0] != negpos || read.negcount > 0) && negpos == fragpos && dir * fragend > dir * (*curr).seen_end) { (*curr).seen_end = fragend; printf("C\n");}
	     //if((dir == -1 || empty(&pvepos) || head(&pvepos).startend[0] != negpos ) &&  negpos == fragpos && read.negcount > 0 && dir * fragend > dir * (*curr).seen_end) { (*curr).seen_end = fragend; printf("C\n");}
	     //if(negpos == fragpos && read.negcount > 0 && dir * fragend > dir * (*curr).seen_end && !read.changechr) { (*curr).seen_end = fragend; printf("C\n");}
	     //if(read.negcount > 0 && dir * fragend > dir * (*curr).seen_end && ((!read.changechr && negpos == fragpos) || (read.changechr && fragpos == (read.pos + par.seqlength - par.fraglength)))) { (*curr).seen_end = fragend; printf("C\n");}
	     if(read.negcount > 0 && dir * fragend > dir * (*curr).seen_end)
	       {
		 if(!read.changechr && negpos == fragpos)
		   { 
		     (*curr).seen_end = fragend; 
		     //printf("C\t%ld\n", (*curr).seen_end);
		   }
		 if(read.changechr && (read.pos + par.seqlength - par.fraglength)) 
		   { 
		     (*curr).seen_end = read.pos + par.seqlength - 1; 
		     //printf("D\t%ld\n", (*curr).seen_end);
		   }
	       }
	   }
	 
	 if(newtot <= 0.00001)  //End of region
	   {
	     //printf("A6aa\n");
	     if(dir == 1)
	       { store_results(read, *possend, par, lastregion);
		 *curr = new_region(read, par, fragpos, fragend, dir); 
		 *possend = *curr;
	       }
	     else{ (*curr).total = newtot; return; }
	   }
	 
	 else
	   {	     		  	      
	     //(*curr).total = newtot;
	     (*curr).total = (par.maxtotal && newtot > par.maxtotal) ? par.maxtotal : newtot;
	     if(newtot >= (*possend).total && read.count > 0)     //New high value of function
	       {
		 //printf("%f %f %ld %d\n", newtot, (*possend).total, read.pos, dir);
		 *possend = *curr;
	       }
	   }
       }
   }


void update_fragpos(long int *fragpos, long int seqlength, long int fraglength, struct readinfo read)
{
  (*fragpos)++;
  long int tempfragpos;
  tempfragpos = read.pos + seqlength - fraglength;
  if(!empty(&pvepos) && head(&pvepos).startend[0] < tempfragpos)
       { tempfragpos = head(&pvepos).startend[0]; }
 
 
  if(read.changechr && *fragpos <= (read.pos + read.changechr - fraglength)  && tempfragpos > (read.pos + read.changechr - fraglength) && read.negcount > 0)
    {
      *fragpos = read.pos + read.changechr - fraglength;
      return;
    }


  if(tempfragpos > *fragpos)
    { *fragpos = tempfragpos; }
}


//***Separately treat negative and positive strand reads***
int separate(char strnd[2], int pairlen, long int pos, int bootstrap)
{
  //printf("S:%s\n",strnd);
  // ***Negative reads ignored at moment for paired end reads***
  if(!strcmp(strnd,negstring)) 
    {  
      return(rpois(bootstrap,0));
      //return(1);
    }
  // ***Positive reads stored in pvepos
  else
    {
    if(!strcmp(strnd,pvestring)) 
      {
	int pois = rpois(bootstrap,0);
	if(pois)
	  { push(pos, pos + pairlen - 1, pois, &pvepos); }
	//push(pos, pos + pairlen - 1, 1, &pvepos);
	return(0);
	//push(pos + pairlength -1, &pveend);
      }
    else{printf("Strand not recognised.\n");exit(1);}
    }
}


struct region new_region(struct readinfo read, struct param par, long int fragpos, long int fragend, int dir) //INCLUDES REF
{
  struct region curr;
  //printf("NR\t%f\n", read.count);
  //if(count > bg + posbg)
  // {
  //printf("NR: %f\n", read.count - par.bg - par.posbg - read.refcount*par.refpen);
  curr.total = (read.count - par.bg - par.posbg - read.refcount * par.refpen > 0) ? (read.count - par.bg - par.posbg - read.refcount * par.refpen) : (0 - par.threshold) ;
  curr.print_total = (read.count - par.bg - par.posbg - read.refcount * par.refpen > 0) ? (read.count - par.bg - par.posbg - read.refcount * par.refpen) : (0 - par.threshold) ;
  if(par.maxtotal && curr.total > par.maxtotal) { curr.total = par.maxtotal; }
  //newtot = total;
  
  strcpy(curr.chr,read.chr);
  curr.start = fragpos;
  curr.seen_start = read.pos;
  //if(dir == 1 && read.pos + read.pairlength - 1 > fragend) {curr.seen_start = read.prevpos; printf("ADJ\n");}
  //if(dir == 1 && fragpos == read.prevpos) {curr.seen_start = read.prevpos; printf("ADJ\n");}
  if(dir == 1 && fragpos == read.prevpos && (read.count > read.negcount || fragpos != (read.pos + par.seqlength - par.fraglength))) {curr.seen_start = read.prevpos;}
  //curr.seen_start = read.pos;
  //curr.hightot = curr.total;
  curr.count = read.count;
  curr.refcount = read.refcount;
  //curr.end = fragpos + dir * (fraglength-1);
  curr.end = fragend;
  //curr.seen_end = read.pos + dir * (read.pairlength - 1);
  curr.seen_end = curr.seen_start + dir * (read.pairlength - 1);
  curr.pos = 1;
  curr.peakheight = 0;
  curr.peakheightmax = 0;
  curr.summit = 0;
  //printf("NRNT: %s\t%ld\t%ld\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%ld\t%ld\t%ld\n", read.chr, read.pos, fragpos, fragend, curr.total, read.count, read.negcount, par.bg, par.posbg, read.refcount, par.refpen, read.prevpos, curr.seen_start, curr.seen_end);
  //num_above_bg=read.count-read.refcount;
      //}

  while(!empty(&peakheights)) { shift(&peakheights); }
  while(!empty(&peaksummits)) { shift(&peaksummits); }

  if(par.ref)
    {
      while(!empty(&refpvepos) && head(&refpvepos).startend[0] < fragpos)
	{ shift(&refpvepos); }
      while(!empty(&refnvepos) && head(&refnvepos).startend[0] < fragpos)
	{ shift(&refnvepos); }
    }
  

  //if(cnt(&refbackpos)){printf("Tail refbackpos: %ld\n", head(&refbackpos).startend[0]);}
  return(curr);
}


 
int set_back_fragpos(long int *fragpos, long int *fragstart)
{
  int backstrand = -1;
  if(!empty(&backpvepos) || !empty(&backnvepos))
    { 
      if(!empty(&backnvepos))
	{
	  *fragpos = head(&backnvepos).startend[1];
	  *fragstart = head(&backnvepos).startend[0];
	  if(!empty(&backpvepos))
	    {
	      if(*fragpos < head(&backpvepos).startend[1])
		{ 
		  *fragpos = head(&backpvepos).startend[1];
		  *fragstart = head(&backpvepos).startend[0];
		  backstrand = 1;
		}
	    }
	}
      else
	{
	  *fragpos = head(&backpvepos).startend[1];
	  *fragstart = head(&backpvepos).startend[0];
	  backstrand = 1;
	}
    }
  else
    { *fragpos = 0;}
  return backstrand;
}


double sample_back_count(long int fragpos, long int *fragstart, double threshold, stack *stk)
{
  double count = 0;
  /*while(!empty(&backpos) && head(&backpos).startend[1] == fragpos)
     {
       if(head(&backpos).startend[0] < *fragstart) { *fragstart = head(&backpos).startend[0]; }
       count += head(&backpos).count;
       push_peakheight(-1, fragpos, head(&backpos).startend[0], head(&backpos).count, back, possend);
       shift(&backpos);
    }
  */
  if(!empty(stk))
    {
      elem *p = stk -> head;
      int i = 0;
     
      if(p -> d[1] == fragpos) 
	{ 
	  count += p -> s; 
	  if(p -> d[0] < *fragstart) { *fragstart = p -> d[0]; }
	  while(i < (stk->cnt - 1)  && p -> next -> d[1] == fragpos)
	    {	      
	      p = p->next;
	      count += p -> s;
	      if(p -> d[0] < *fragstart) { *fragstart = p -> d[0]; }
	      i++;
	    }
	}      
      // printf("SBC\t%d\t%f\n", i, count);
    }
  
  if(count > threshold){ count = threshold;}
  return(count);
}


void output_regions(struct region possend, struct region *olr, struct region pvestrandregion, struct param par)
{
  int tempstart = possend.end;
  possend.end = possend.start;
  possend.start = tempstart;
  //printf("ORS:%f\t%f\t%f\n", possend.count, possend.refcount, possend.total);
  if((possend.count - possend.refcount*par.refpen ) >= par.min_above_bg && possend.total >= par.resultcutoff)
    {
      calculate_peak_summit(&possend, head(&peakheights).startend[0], tail(&peakheights).startend[0], &possend);

      /*
	int i;
	for(i = 0; i < cnt(&peaksummits); i++)
	 {
	printf("%ld ", (long int)return_position(i, &peaksummits));
	 }
	printf("\n");
      */
      
      possend.summit = -1 * median(&peaksummits);
      // fprintf(out_fp,"%s\t%ld\t%ld\t%f\t%ld\t%d\t%f\t%f\t%ld\t%ld\t-\n", possend.chr, possend.start, possend.end, possend.count, possend.end-possend.start+1, possend.pos, possend.total,possend.refcount, possend.seen_end, possend.seen_start);
      //printf("%s\t%ld\t%ld\t%f\t%ld\t%d\t%f\t%f\t%f\t%f\t%ld\t%ld\t-\n", possend.chr, possend.start, possend.end, possend.count, possend.end-possend.start+1, possend.pos, possend.total, possend.refcount, possend.peakheightmax, possend.summit, possend.seen_end, possend.seen_start);
      // while(!empty(&peaksummits)) { printf("%f ", shift(&peaksummits).count); } printf("\n");      );
      //printf("OR\n");

      
      long int temp_possend_seenend = possend.seen_end;
      possend.seen_end = possend.seen_start;
      possend.seen_start = temp_possend_seenend;
      if(!(pvestrandregion.total) || (possend.start > pvestrandregion.end))
	{
	  //printf("O1\t%f\t%ld\t%ld\n",pvestrandregion.total, possend.start, pvestrandregion.end);
	  //printf("O1a\t%f\t%ld\t%ld\t%ld\t%ld\n",(*olr).total, (*olr).start, (*olr).end, (*olr).seen_start, (*olr).seen_end);
	  unshiftregion(possend, &outregions);
	}
      else
	{
	  //printf("O1b\t%f\t%ld\t%ld\t%ld\t%ld\n",(*olr).total, (*olr).start, (*olr).end, (*olr).seen_start, (*olr).seen_end);
	  if((possend.start <= pvestrandregion.end) && (possend.end >= pvestrandregion.start))
	    {
	      //printf("O2\n");
	      if((*olr).start > possend.start) { (*olr).start = possend.start; }
	      if((*olr).end < possend.end) { (*olr).end = possend.end; }	
	      if((*olr).seen_start > possend.seen_start) { (*olr).seen_start = possend.seen_start; }
	      if((*olr).seen_end < possend.seen_end) { (*olr).seen_end = possend.seen_end; }  
	      if((*olr).count < possend.count) { (*olr).count = possend.count; }
	      //(*olr).length = (*olr).end - (*olr).start + 1;
	      if((*olr).pos < possend.pos) { (*olr).pos = possend.pos; }
	      if((*olr).refcount < possend.refcount) { (*olr).refcount = possend.refcount; }
	      if((*olr).total < possend.total) { (*olr).total = possend.total; }
	      if((*olr).print_total < possend.print_total) { (*olr).print_total = possend.print_total; }
	      if((*olr).peakheightmax < possend.peakheightmax) 
		{ 
		  (*olr).peakheightmax = possend.peakheightmax; 
		  (*olr).summit = possend.summit;
		}	      
	      //printf("O2a\t%f\t%ld\t%ld\n",(*olr).total, (*olr).start, (*olr).end);
	    }
	  else
	    {
	      if((*olr).total)
		{
		  //printf("O3\n");	       
	      unshiftregion(*olr, &outregions);
	      (*olr).total = 0;
	      unshiftregion(possend, &outregions);
		}
	      else
		{
		  //printf("O4\n");
		  unshiftregion(possend, &outregions);
		}
	    }
	}
    } 
  //return(olr);
}


struct region check_and_print_regions(struct region olr, struct region *lastregion, int comp, double binom_p, int peakversion)
  {
    if(olr.total)
      { 
	//printf("OH\n");
	unshiftregion(olr, &outregions);
      }



    while(!emptyregion(&outregions))
      { 
	struct region print_region = shiftregion(&outregions);
	//printf("%s\t%ld\t%ld\t%ld\t%d\t%f\t%f\t%ld\t%ld\t%ld\tK\n", print_region.chr, print_region.start, print_region.end, print_region.end-print_region.start+1, print_region.pos, print_region.print_total, print_region.peakheightmax, (long int)floor(print_region.summit), print_region.seen_start, print_region.seen_end);
	if(strcmp(print_region.chr,(*lastregion).chr) || (*lastregion).end < print_region.start)
	  {
	    if((*lastregion).total)
	      { 
		if(comp) { compare_lastregion(lastregion);}
		else {print_out_lastregion(lastregion, binom_p, peakversion); }
	      }
	    *lastregion = print_region;	  
	  }
	else
	  {
	    (*lastregion).end = print_region.end;
	    (*lastregion).count += print_region.count;
	    (*lastregion).pos += print_region.pos;
	    (*lastregion).refcount += print_region.refcount;
	    (*lastregion).seen_end = print_region.seen_end;
	    if(print_region.total > (*lastregion).total) { (*lastregion).total = print_region.total; }
	    if(print_region.print_total > (*lastregion).print_total) { (*lastregion).print_total = print_region.print_total; }
	    if(print_region.peakheightmax > (*lastregion).peakheightmax) 
	      { 
		(*lastregion).peakheightmax = print_region.peakheightmax;  
		(*lastregion).summit = print_region.summit;
	      }
	  }

      }
    
    //printf("%s\t%ld\t%ld\t%ld\t%d\t%f\t%f\t%ld\n", (*lastregion).chr, (*lastregion).start, (*lastregion).end, (*lastregion).end-(*lastregion).start+1, (*lastregion).pos, (*lastregion).print_total, (*lastregion).peakheightmax, (long int)floor((*lastregion).summit));
    return (*lastregion);
  }


void poiscalc(double ratio, double ppois[10], double refppois[10])
{
  double refppoi[] = { 0.3678794, 0.7357589, 0.9196986, 0.9810118, 0.9963402, 
  	       0.9994058, 0.9999168, 0.9999898, 0.9999989, 0.9999999};
  int i;
  for(i =0; i < 10; i++) {refppois[i] = refppoi[i];}

  //double ppois[10];
  double fac[] = {1,1,2,6,24,120,720,5040,40320,362880}; 

  ppois[0] = exp(ratio * -1);
  for(i = 1; i < 10; i++)
    {
      ppois[i] = ppois[i-1] + pow(ratio,i)*exp(ratio * -1)/fac[i];
    }
  return;
}


int rpois(int bootstrap, int ref)
{
  int count = 1;
  if(bootstrap)
    {
      double p = (double)rand() / ((double)(RAND_MAX)+(double)(1));
      count = 0;
      if(ref)
	{
	  while(p > refppois[count] && count++ < 10)
	    {}
	}
      else
	{
	  while(p > ppois[count] && count++ < 10)
	    {}
	}
    }
  //printf("%d\n", count);
  return count;
}

double pbinom(int q, int n, double p)
{
  double t = ((double)q - (double)n*p)/sqrt((double)n*p*(1-p));
  //printf("%d\t%d\t%f\t%f\t%f\n", q, n, p, t, pnorm(t));
  return pnorm(t);
}

double pnorm(double x)
{
  return 0.5*(1+erf(x/sqrt(2)));
}
