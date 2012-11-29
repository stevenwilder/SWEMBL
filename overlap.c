#include "SWEMBL.h"

void compare_lastregion(struct region *lastregion)
{
  overlap.nseqregions++;
  overlap.lenseqregions += (*lastregion).end - (*lastregion).start; // + 1;
  //printf("%s\t%ld\t%ld\n", (*lastregion).chr,  (*lastregion).start, (*lastregion).end);

  int compendfile = 0;
  //printfshiftchr("N:%s\t%s\t%ld\t%ld\n", lastchrelem(&compchrs)->chr, (*lastregion).chr, emptychr(&compchrs) ? 0 : headchr(&compchrs).startend[1], (*lastregion).start);
  char prevcompchr[50];
  strcpy(prevcompchr,lastchrelem(&compchrs)->chr);
  //printf("LC:%s\n", prevcompchr);
  
  compendfile = store_comp_location(prevcompchr);
  

  while((emptychr(&compchrs) || strcmp(firstchrelem(&compchrs)->chr, (*lastregion).chr) || empty(firstchr(&compchrs)) || headchr(&compchrs).startend[1] <= (*lastregion).start) && !(compendfile && emptychr(&compchrs)))
    {
      
      
      int seen = seenchr((*lastregion).chr, &compchrs);
      // printf("W:%s\t%s\t%ld\t%ld\t%d\n", lastchrelem(&compchrs)->chr,firstchrelem(&compchrs)->chr, (*lastregion).chr,  (*lastregion).start, headchr(&compchrs).startend[0],seen);
      // printf("W:%s\t%s\t%s\t%ld\n", lastchrelem(&compchrs)->chr,firstchrelem(&compchrs)->chr, (*lastregion).chr,  (*lastregion).start);
      //printf("W:%s\t%s\t%s\t%ld\t%ld\t%d\t%d\n", lastchrelem(&compchrs)->chr, firstchrelem(&compchrs)->chr, (*lastregion).chr, (emptychr(&compchrs) || empty(firstchr(&compchrs))) ? 0 : headchr(&compchrs).startend[1], (*lastregion).start, compendfile, seen);

      //  ##Comp regions are on earlier chromosomes
      while(!emptychr(&compchrs) && strcmp(firstchrelem(&compchrs)->chr, (*lastregion).chr) && seen)
	{
	  //printf("RR\n");
	  while(!empty(firstchr(&compchrs)))
	    {
	      //printf("s1\n");
	      if(overlap.print) {nearestregion(&compchrs, lastregion);}
	      shiftchr(&compchrs);
	    }
	  if(!emptychr(&compchrs)){removechr(&compchrs);}	  
	  overlap.compprioroverlap = 0;
	}

      //  ##Comp regions are earlier in chromosome
      while(!emptychr(&compchrs) && !empty(firstchr(&compchrs)) && !(!strcmp(firstchrelem(&compchrs)->chr, (*lastregion).chr) && headchr(&compchrs).startend[1] > (*lastregion).start)  && seen) 
	{
	  //printf("s2\n");
	  if(overlap.print && !overlap.compprioroverlap){nearestregion(&compchrs, lastregion);}
	  shiftchr(&compchrs);	  
	  overlap.compprioroverlap = 0;
	}
      //printf("ASC\n");
      
      compendfile = store_comp_location(prevcompchr);
      if(compendfile){break;}
      //removechr(&compchrs);
    }
  
  /*   while(!emptychr(&compchrs) && !strcmp(firstchrelem(&compchrs)->chr, (*lastregion).chr) && !empty(firstchr(&compchrs)) && (*lastregion).start <= headchr(&compchrs).startend[1] && headchr(&compchrs).startend[0] <= (*lastregion).end)
    {
    printf("%s\t%s\t%ld\t%ld\t%ld\t%ld\n", firstchrelem(&compchrs)->chr, (*lastregion).chr, headchr(&compchrs).startend[0], headchr(&compchrs).startend[1], (*lastregion).start, (*lastregion).end);
      
      if(!empty(firstchr(&compchrs))) {shiftchr(&compchrs);}
      store_comp_location(prevcompchr);
      } */
  
  int lastregionoverlap = 0;
  while(!emptychr(&compchrs) && !strcmp(firstchrelem(&compchrs)->chr, (*lastregion).chr) && !empty(firstchr(&compchrs)) && headchr(&compchrs).startend[0] <= (*lastregion).end)
    {
      lastregionoverlap = 1;
      overlap.ncompregionsoverlap += 1 - overlap.compprioroverlap;
      overlap.compprioroverlap = 0;

      overlap.lenoverlap += 1;
      overlap.lenoverlap += (headchr(&compchrs).startend[1] < (*lastregion).end) ? 
	                         headchr(&compchrs).startend[1] : (*lastregion).end;
      overlap.lenoverlap -= (headchr(&compchrs).startend[0] > (*lastregion).start) ?
	                         headchr(&compchrs).startend[0] : (*lastregion).start;
      //printf("LO:%ld%s\t%ld\t%ld\t%f\t%ld\t%d\t%f\t%f\t%s\t%ld\t%ld\n", overlap.lenoverlap, (*lastregion).chr, (*lastregion).start, (*lastregion).end, (*lastregion).count, (*lastregion).end-(*lastregion).start+1, (*lastregion).pos, (*lastregion).refcount,(*lastregion).total, firstchrelem(&compchrs)->chr, headchr(&compchrs).startend[0], headchr(&compchrs).startend[1]);
      if(overlap.print)
	{
	  fprintf(overreg_fp, "%s\t%ld\t%ld\t", (*lastregion).chr, (*lastregion).start, (*lastregion).end+1);
	  if(!overlap.mode){ fprintf(overreg_fp, "%f\t%ld\t%d\t%f\t%f\t", (*lastregion).count, (*lastregion).end-(*lastregion).start, (*lastregion).pos, (*lastregion).refcount,(*lastregion).print_total); }
	  fprintf(overreg_fp, "%s\t%ld\t%ld\t0\n", firstchrelem(&compchrs)->chr, headchr(&compchrs).startend[0], headchr(&compchrs).startend[1]);
	}
      if((*lastregion).end <= headchr(&compchrs).startend[1]) 
	{ overlap.compprioroverlap = 1;}
      if(!empty(firstchr(&compchrs))) 
	{
	  //printf("s3\n");
	  if(overlap.print)
	    {
	      fprintf(overcomp_fp, "%s\t%ld\t%ld\t%s\t%ld\t%ld\t", firstchrelem(&compchrs)->chr, headchr(&compchrs).startend[0], headchr(&compchrs).startend[1], (*lastregion).chr, (*lastregion).start, (*lastregion).end+1);
	      if(!overlap.mode){ fprintf(overcomp_fp, "%f\t%ld\t%d\t%f\t%f\t", (*lastregion).count, (*lastregion).end-(*lastregion).start, (*lastregion).pos, (*lastregion).refcount,(*lastregion).print_total); }
	      fprintf(overcomp_fp, "0\n");
		
	    }
	  store_comp_location(prevcompchr);

	  if(cnt(firstchr(&compchrs)) > 1 && nextheadchr(&compchrs).startend[0] <= (*lastregion).end)
	    {	    
	      shiftchr(&compchrs);
	      overlap.compprioroverlap = 0;
	    }
	  else
	    {break;}
	}
      //store_comp_location(prevcompchr);
    }
  overlap.nseqregionsoverlap += lastregionoverlap;
  if(overlap.print && !lastregionoverlap)
    {
      //printf("%s\t%ld\t%ld\t%s\t%ld\t%ld\n",(*lastregion).chr, (*lastregion).start, (*lastregion).end, firstchrelem(&compchrs)->chr, headchr(&compchrs).startend[0], headchr(&compchrs).startend[1] );
      if(!emptychr(&compchrs) && !strcmp(firstchrelem(&compchrs)->chr, (*lastregion).chr) && !empty(firstchr(&compchrs)) &&(*lastregion).end <= headchr(&compchrs).startend[0])
	{
	  if(strcmp(prev.compchr, (*lastregion).chr) || (headchr(&compchrs).startend[0] - (*lastregion).end) <= ((*lastregion).start - prev.compend))
	    { 
	      fprintf(overreg_fp, "%s\t%ld\t%ld\t", (*lastregion).chr, (*lastregion).start, (*lastregion).end+1);
	      if(!overlap.mode){ fprintf(overreg_fp, "%f\t%ld\t%d\t%f\t%f\t", (*lastregion).count, (*lastregion).end-(*lastregion).start, (*lastregion).pos, (*lastregion).refcount,(*lastregion).print_total); }
	      fprintf(overreg_fp, "%s\t%ld\t%ld\t%ld\n", firstchrelem(&compchrs)->chr, headchr(&compchrs).startend[0], headchr(&compchrs).startend[1], headchr(&compchrs).startend[0] - (*lastregion).end);
	    }
	  else
	    {
	      fprintf(overreg_fp, "%s\t%ld\t%ld\t", (*lastregion).chr, (*lastregion).start, (*lastregion).end+1);
	      if(!overlap.mode){ fprintf(overreg_fp, "%f\t%ld\t%d\t%f\t%f\t", (*lastregion).count, (*lastregion).end-(*lastregion).start, (*lastregion).pos, (*lastregion).refcount,(*lastregion).print_total); }
	      fprintf(overreg_fp, "%s\t%ld\t%ld\t%ld\n", prev.compchr, prev.compstart, prev.compend, (*lastregion).start - prev.compend);
	    }
	}
      else
	{
	  if(!strcmp(prev.compchr, (*lastregion).chr))
	    {
	      fprintf(overreg_fp, "%s\t%ld\t%ld\t", (*lastregion).chr, (*lastregion).start, (*lastregion).end+1);
	      if(!overlap.mode){ fprintf(overreg_fp, "%f\t%ld\t%d\t%f\t%f\t", (*lastregion).count, (*lastregion).end-(*lastregion).start, (*lastregion).pos, (*lastregion).refcount,(*lastregion).print_total); }
	      fprintf(overreg_fp, "%s\t%ld\t%ld\t%ld\n", prev.compchr, prev.compstart, prev.compend, (*lastregion).start - prev.compend);
	    }
	  else
	    { 
	      fprintf(overreg_fp, "%s\t%ld\t%ld\t", (*lastregion).chr, (*lastregion).start, (*lastregion).end+1);
	      if(!overlap.mode){ fprintf(overreg_fp, "%f\t%ld\t%d\t%f\t%f\t", (*lastregion).count, (*lastregion).end-(*lastregion).start, (*lastregion).pos, (*lastregion).refcount,(*lastregion).print_total); }
	      fprintf(overreg_fp, "NF\tNA\tNA\tNA\n");
	    }
	}
    }
  prev.last_but_1_region = *lastregion;
  //printf("CLRE\n");
}



int store_comp_location(char prevcompchr[50])
{
  char compchr[50];
  //int compn, compchrnum = -1;
  char compline[1000];

  if(fgets (compline, 100, comp_fp) != NULL)
    {
      //printf("%s", compline);
      int l = 0;
      char *compsplit = NULL;
      
      compsplit = strtok( compline, "\t" );
      long int compstart = 0;
      //printf("%s\n", compsplit);
      
      while(l < 3)
	{
	  switch(l)
	    {
	    case(0): 
	      {
		strcpy(compchr, compsplit); 
		if(strcmp(compchr,prevcompchr))
		{
		  addchr(compchr, &compchrs);
		}
	      strcpy(prevcompchr,compchr);
	      //printf("PC:%s\t%s\n",prevcompchr,compchr);
	      break;
	      }
	      //default: { compresults[compchrnum][compn][l-1] = atoi(compsplit); break; }
	    case(1): compstart = atoi(compsplit); break;
	    case(2): 
	      {
		overlap.lencompregions += atoi(compsplit) - compstart;// + 1;
		//push(compstart, atoi(compsplit), compresults[compchrnum]); break;
		pushchr(compstart, atoi(compsplit), 1, &compchrs);
	      }
	    }
	  l++;	  
	  compsplit = strtok( NULL, "\t" );	  	  
	  //printf("%s\n", compline);
	  //if(compsplit != NULL){printf("%s\n", compsplit);}
	} 
      // printf("%d\t%d\t%ld\t%ld\n",compchrnum, compn, tail(compresults[compchrnum]).startend[0], tail(compresults[compchrnum]).startend[1]);
      overlap.ncompregions++;
      return(0);
    }
  else
    { return(1); }
  //overlap.nchr = compchrnum+1;
}


void nearestregion(struct chrstack *cstk, struct region *lastregion)
{
  //printf("NER:%s\t%ld\t%ld\t%s\t%ld\t%ld\t%s\t%ld\t%ld\n", firstchrelem(&compchrs)->chr, headchr(&compchrs).startend[0], headchr(&compchrs).startend[1], (*lastregion).chr, (*lastregion).start, (*lastregion).end, prev.last_but_1_region.chr, prev.last_but_1_region.start, prev.last_but_1_region.end);

  if(!strcmp(prev.last_but_1_region.chr, firstchrelem(&compchrs)->chr) && prev.last_but_1_region.end > headchr(&compchrs).startend[0])
    {
      //printf("MISS: %s\t%ld\t%s\t%ld\n", firstchrelem(&compchrs)->chr, headchr(&compchrs).startend[0], prev.last_but_1_region.chr, prev.last_but_1_region.end);
      return;
    }

  if(!strcmp(firstchrelem(cstk)->chr, (*lastregion).chr) && headchr(&compchrs).startend[1] < (*lastregion).start)
    {
      if(strcmp(prev.last_but_1_region.chr, firstchrelem(cstk)->chr) || ((*lastregion).start - headchr(&compchrs).startend[1]) <= (headchr(&compchrs).startend[0] - prev.last_but_1_region.end))
	{
	  fprintf(overcomp_fp, "%s\t%ld\t%ld\t%s\t%ld\t%ld\t", firstchrelem(&compchrs)->chr, headchr(&compchrs).startend[0], headchr(&compchrs).startend[1], (*lastregion).chr, (*lastregion).start, (*lastregion).end+1);
	  if(!overlap.mode){ fprintf(overcomp_fp, "%f\t%ld\t%d\t%f\t%f\t", (*lastregion).count, (*lastregion).end-(*lastregion).start, (*lastregion).pos, (*lastregion).refcount,(*lastregion).print_total); }
	  fprintf(overcomp_fp, "%ld\n", (*lastregion).start - headchr(&compchrs).startend[1]);
	}
      else
	{
	  fprintf(overcomp_fp, "%s\t%ld\t%ld\t%s\t%ld\t%ld\t", firstchrelem(&compchrs)->chr, headchr(&compchrs).startend[0], headchr(&compchrs).startend[1], prev.last_but_1_region.chr, prev.last_but_1_region.start, prev.last_but_1_region.end+1);
	  if(!overlap.mode){ fprintf(overcomp_fp, "%f\t%ld\t%d\t%f\t%f\t", prev.last_but_1_region.count, prev.last_but_1_region.end-prev.last_but_1_region.start, prev.last_but_1_region.pos, prev.last_but_1_region.refcount, prev.last_but_1_region.print_total); }
	  fprintf(overcomp_fp, "%ld\n", headchr(&compchrs).startend[0] - prev.last_but_1_region.end);
	}
    }
  
  else
    {
      if(!strcmp(prev.last_but_1_region.chr, firstchrelem(cstk)->chr))
	{
	  fprintf(overcomp_fp, "%s\t%ld\t%ld\t%s\t%ld\t%ld\t", firstchrelem(&compchrs)->chr, headchr(&compchrs).startend[0], headchr(&compchrs).startend[1], prev.last_but_1_region.chr, prev.last_but_1_region.start, prev.last_but_1_region.end+1);
	  if(!overlap.mode){ fprintf(overcomp_fp, "%f\t%ld\t%d\t%f\t%f\t", prev.last_but_1_region.count, prev.last_but_1_region.end-prev.last_but_1_region.start, prev.last_but_1_region.pos, prev.last_but_1_region.refcount, prev.last_but_1_region.print_total); }
	  fprintf(overcomp_fp, "%ld\n", headchr(&compchrs).startend[0] - prev.last_but_1_region.end);
	}
      else
	{
	  fprintf(overcomp_fp, "%s\t%ld\t%ld\tNF\tNA\tNA\tNA\n", firstchrelem(&compchrs)->chr, headchr(&compchrs).startend[0], headchr(&compchrs).startend[1]);
	}
    }
}



void print_out_overlap()
{
  fprintf(out_fp, "%d\t%ld\t%d\t%ld\t%d\t%ld\t%d\n", overlap.nseqregions, overlap.lenseqregions, overlap.nseqregionsoverlap, overlap.lenoverlap, overlap.ncompregions, overlap.lencompregions, overlap.ncompregionsoverlap);
}
