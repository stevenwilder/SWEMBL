#include "SWEMBL.h"

void establish_wiggle(char infile[1000], char contigfile[1000], struct contigconv conversion[], int *wigonevalue)
{
  fprintf(wig_fp,"track name=\"%s\" description=\"peaks %s\" type=\"wiggle_0\" color=50,50,150 yLineMark=0.0 yLineOnOff=off visibility=2 priority=1 maxHeightPixels=40:40:2\n", infile, infile);
  *wigonevalue = 0;
  

  if(strcmp(contigfile,""))
    {    
      int c=0;
      char contigline[1000];
      
      while(egets(contigline, 1000, NT_fp) != NULL)
	{
	  char *nt;
	  nt = strchr(contigline, '\n');
	  if(nt){*nt = '\0';}
	
	  char *ntresult = NULL;
	  ntresult = strtok(contigline, "\t");
	  int i = 0;
	    
	  while( ntresult != NULL)
	    {
	      switch(i)
		{
		case(0): {strcpy(conversion[c].contig, ntresult); break;}
		case(1): { 
		  char nt3[] = {ntresult[0], ntresult[1], ntresult[2], '\0'};
		//printf("NT3:%s\n", nt3);
		  if(!strcmp(nt3,"chr"))
		    {
		      int c;
		      char rmchr[strlen(ntresult)-3];
		      for(c=0; c<= (strlen(ntresult) -3);c++)
			{
			  rmchr[c] = ntresult[c+3];
			}
		      //printf("rmchr:%s\n", rmchr);
		      strcpy(ntresult, rmchr);
		    }
		  strcpy(conversion[c].chr_random, ntresult); 
		  break;}
		case(2): {conversion[c].start = atoi(ntresult); break;}
		case(3): {conversion[c++].end = atoi(ntresult); break;}
		}
	      i++;
	      ntresult = strtok(NULL, "\t");
	    }
	}
    }
}


void push_wigcounts(struct readinfo read, long int fragpos, long int prevfragpos, long int fragend, int fraglength, struct contigconv conversion[], int contigrows, char wigchr[], long int *wigstart, int paired, int *wigonevalue)
{
  //if(strcmp(read.chr,"10"))
    {
      //printf("PW:%ld\t%ld\t%ld\t%s\t%s\t%s\t%ld\n", fragpos, fragend, prevfragpos, read.chr, read.prevchr, wigchr, *wigstart);
    }
  if(!empty(&wigcounts))
    {
      int w;
      /*if(!paired)
	{
	  if(strcmp(read.chr,read.prevchr) || fragpos - prevfragpos > fraglength)
	    { export_wiggle(conversion, contigrows, wigchr, wigstart, fraglength, paired);}
	  else
	    {
	      for(w=0; w<(fragpos-prevfragpos-1); w++)
		{ push(0, 0, &wigcounts); }
	      push(1, 0, &wigcounts);
	    }
	}
	else*/
	{
	  if(!(!strcmp(wigchr,read.chr) && fragpos == *wigstart) && (strcmp(read.chr,wigchr) || fragpos >= *wigstart + cnt(&wigcounts)))
	    {export_wiggle(conversion, contigrows, wigchr, wigstart, fraglength, paired, *wigonevalue);}
	  else
	    { 
	      //printf("sdf\n");
	      /*for(w=*wigstart + cnt(&wigcounts); w <= fragend+1; w++)
		{ push(0, 0, 1, &wi gcounts);}
	      add_to_position(1, 0,  fragpos-*wigstart, &wigcounts);    ***UPDATED for double
	      add_to_position(1, 1, fragend-*wigstart+1, &wigcounts);*/

	      for(w=*wigstart + cnt(&wigcounts); w <= fragend+1; w++)
		{ push(0, 0, 0, &wigcounts);}
	      add_to_position(read.count,  fragpos-*wigstart, &wigcounts);
	      add_to_position(-1*read.count, fragend-*wigstart+1, &wigcounts);
	      *wigonevalue = 1;
	    }
	}

      /*      if(!(!strcmp(wigchr,read.chr) && fragpos == *wigstart) && (strcmp(read.chr,read.prevchr) || fragpos - prevfragpos > fraglength))
	      { export_wiggle(conversion, contigrows, wigchr, wigstart, fraglength, paired);}
	      else
	      {//Multiple paired end starting same pos
	      int w;
	      if(!paired)
	      {
	      for(w=0; w<(fragpos-prevfragpos-1); w++)
	      { push(0, 0, &wigcounts); }
	      push(1, 0, &wigcounts);
	      }
	      if(paired)
	      { 
	      //for(w=fragpos-prevfragpos+1; w<fraglength; w++)
	       for(w=*wigstart + cnt(&wigcounts); w <= fragend+1; w++)
		 { push(0, 0, &wigcounts);}
	       add_to_position(1, 0, fragpos-*wigstart, &wigcounts);
	       add_to_position(1, 1, fragend-*wigstart+1, &wigcounts);
	     }		 
	     } */
    }
  
  if(empty(&wigcounts))
    {
      *wigonevalue = 0;
      strcpy(wigchr, read.chr);
      if(fragpos > 1)
	{ *wigstart = fragpos;}
      else { *wigstart = 1; }
      //push(1,0, 1, &wigcounts);
      push(0, 0, read.count, &wigcounts);
      //if(paired)
	{ 
	  int w;
	  for(w=*wigstart+1; w<=fragend; w++)
	    { push(0, 0, 0, &wigcounts); }
	    //{ push(0, 0, 1, &wigcounts);}
	  //push(0, 1, 1, &wigcounts);
	  push(0, 0, -1*read.count, &wigcounts);
	  //add_to_position(1, 1, fragend-*wigstart, &wigcounts);
	}		 
    }
}


void export_wiggle(contigconv conv[], int n, char wigchr[], long int *wigstart, int fraglength, int paired, int wigonevalue)
{
  //NT
  /*
  int r;
  for(r = 0; r < n ; r++)
    {
      if(!strcmp(wigchr, conv[r].contig))
	{
	  strcpy(wigchr, conv[r].chr_random);
	  *wigstart += conv[4].start - 1;
	  *wigstart = (*wigstart < 1) ? 1 : *wigstart;
	  break;
	}
    }
  */
  
  
  //printf("EW\n");
  if(!strcmp(wigchr,"MT"))
    {
      strcpy(wigchr,"M");
    }



  int wiglength = cnt(&wigcounts);
  /*if(!paired)
    {
      if(wiglength == 1)
	{
	  fprintf(wig_fp, "fixedStep chrom=chr%s start=%ld step=%d\n%ld\n", wigchr, *wigstart, fraglength, head(&wigcounts).startend[0]);
	}
      
      else
	{
	  fprintf(wig_fp, "fixedStep chrom=chr%s start=%ld step=1\n", wigchr, *wigstart);
	  int r;
	  for(r=0; r<fraglength-1; r++)
	    {
	      push(0,0, &wigcounts);
	    }
	  
	  int sum=0;
	  stack wigarray;
	  initialize(&wigarray);
	  
	  
	  for(r=0; r <= wiglength + fraglength - 2; r++)
	    {
	      if(!paired)
		{
		  if(r < wiglength)
		    {
		      sum += head(&wigcounts).startend[0];
		      push(head(&wigcounts).startend[0], shift(&wigcounts).startend[1], &wigarray);
		    }	      
		  sum -= ((r >= fraglength) ? shift(&wigarray).startend[0] : 0);
		  fprintf(wig_fp,"%d\n",sum);
		}
	    }
	  while(!empty(&wigarray))
	    { shift(&wigarray); }
	}
      while(!empty(&wigcounts))
	{ shift(&wigcounts); }
    }
    else*/
  //{
    // printf("WPOV:%d\n",wiggle_paired_one_value(&wigcounts));


  //if(wiggle_paired_one_value(&wigcounts))
  if(wigonevalue == 0)
    { 
      if(strcmp(filetype,"C"))
	{
	  int wigcnt = floor(head(&wigcounts).count + 0.5);
	  fprintf(wig_fp, "fixedStep chrom=chr%s start=%ld step=%d\n%d\n", wigchr, *wigstart, cnt(&wigcounts)-1, wigcnt);
	}
      else
	{
	  fprintf(wig_fp, "fixedStep chrom=chr%s start=%ld step=%d\n%f\n", wigchr, *wigstart, cnt(&wigcounts)-1, shift(&wigcounts).count);
	}
    }
      
	     
      //fprintf(wig_fp, "fixedStep chrom=chr%s start=%ld step=%d\n%ld\n", wigchr, *wigstart, cnt(&wigcounts)-1, head(&wigcounts).startend[0]);
      ////printf("fixedStep chrom=chr%s start=%ld step=%d\n%ld\n", wigchr, *wigstart, cnt(&wigcounts)-1, head(&wigcounts).startend[0]);
  
  else
    {
      fprintf(wig_fp, "fixedStep chrom=chr%s start=%ld step=1\n", wigchr, *wigstart);
      //printf("fixedStep chrom=chr%s start=%ld step=1\n", wigchr, *wigstart);
      int r; 
      if(strcmp(filetype,"C"))
	{
	  int sum = 0;
	  for(r = 0; r < wiglength - 1; r++)
	    {
	      sum += floor(shift(&wigcounts).count+0.5);
	      fprintf(wig_fp,"%d\n",sum);
	    }
	}
      
      else
	{      	      
	  double sum = 0;
	  int r;
	  for(r = 0; r < wiglength - 1; r++)
	    {
	      //sum += head(&wigcounts).startend[0];
	      //sum -= shift(&wigcounts).startend[1];
	      sum += shift(&wigcounts).count;
	      fprintf(wig_fp,"%f\n",sum);
	      //printf("%f\n",sum);
	    }
	}    
    }
  
  
  while(!empty(&wigcounts))
    { shift(&wigcounts); }
}
   

