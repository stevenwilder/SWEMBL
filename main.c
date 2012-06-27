#include "SWEMBL.h"

int main(int argc, char **argv)
{
  struct param par;
  struct region curr, possend, lastregion;
  lastregion.total = 0; strcpy(lastregion.chr, ""); 
  lastregion.print_total = 0;
  char infile[1000] = "";
  char outfile[1000] = "";
  char reffile[1000] = "";
  char wigfile[1000] = "";
  char contigfile[1000] = "";
  char compfile[1000] = "";
  char compoverlapfile[1000] = "";
  char regoverlapfile[1000] = "";
  int endfile = 0;
  int refendfile = 0;
  char nextline[1000];
  struct readinfo read, refread;
  long int fragpos, prevfragpos = 0;
  long int fragend;
  long int negpos = 0;
  long int refnegpos = 0; 
  long int wigstart;
  int wigonevalue = 0;
  char wigchr[50] = "";
  strcpy(refread.chr,"");
  strcpy(read.prevchr,""); read.count = 0; read.refcount = 0;
  read.pairflag = 18; refread.pairflag = 18; read.qual = 80;
  read.changechr = 0; refread.qual = 80;
  //strcpy(filetype,"M");
  strcpy(filemode, "w");
  strcpy(pvestring, "+");
  strcpy(negstring, "-");
  strcpy(read.strand, pvestring);
  
  overlap.print = 0;
  overlap.mode = 0;
  overlap.nseqregions = 0;
  overlap.nseqregionsoverlap = 0;
  overlap.ncompregions = 0;
  overlap.ncompregionsoverlap = 0;
  overlap.lenseqregions = 0;
  overlap.lenoverlap = 0;
  overlap.lencompregions = 0;
  overlap.compprioroverlap = 0;
  rseed = 0;

  initialize(&pvepos);
  initialize(&backpvepos);
  initialize(&backnvepos);
  initializeregion(&outregions);
  initialize(&refpvepos);
  initialize(&refnvepos);
  initialize(&savedpverefs);
  initialize(&savednegrefs);
  initialize(&savedrefpos);
  initialize(&refbackpos);
  initialize(&wigcounts);
  initializechr(&refchrs_pve);
  initializechr(&refchrs_nve);
  initializechr(&compchrs);
  initialize(&peakheights);
  initialize(&peaksummits);



  par = get_parameters(infile, outfile, reffile, wigfile, contigfile, compfile, compoverlapfile, regoverlapfile, argc,argv);
  //printf("GP\n");
  open_files(par.ref, par.wig, par.comp, infile, outfile, reffile, wigfile, contigfile, compfile, compoverlapfile, regoverlapfile);
  //printf("FT:%s\n", filetype);
  endcol = platform_strand_endline();

  if(rseed > 0) { srand(rseed); } else{ srand(time(NULL)); }
  (par.ref && par.totalreads && par.totalrefreads) ? poiscalc((double)par.totalreads / (double)par.totalrefreads, ppois, refppois) : poiscalc(1, ppois, refppois);

  if(par.ref && par.totalreads && par.totalrefreads) { par.binom_p = (double)par.totalreads / (double)(par.totalreads + par.totalrefreads); }
  //printf ("p: %f\n", par.binom_p);

  par = get_firstline(&endfile, nextline, par);

  printf("%d total reads.\n", par.totalreads); 
  if(par.ref)
    { printf("%d total reference reads.\n", par.totalrefreads); }

  if(!par.comp)
    { print_output_header(par, infile, reffile); }
  
  int nsavedrefchrs = 0;
  ///struct savedreflengths   savedrefs[par.maxsavedrefchrs];
  struct contigconv conversion[par.contigrows];
  //int totalreadsleft = par.totalreads;
  int totalreadsleft = par.nlines;
  int nreadsleft = floor(par.proportion * totalreadsleft + 0.5);
  //srand(time(NULL));;
  //srand(1);
  
  unsigned int seqlength;
  int fraglength;
  
  strcpy(prev.last_but_1_region.chr,"");
  prev.last_but_1_region.total = 0;
  prev.last_but_1_region.print_total = 0;
  strcpy(prev.compchr,"");

  if(par.wig){ establish_wiggle(infile, contigfile, conversion, &wigonevalue);}

  if(par.ref){ readrefline(&refread, &refendfile, par.reffraglength, par.paired, par.bootstrap); } //PAired ends
  
  if(par.ref)
    {
      //printf("RRF:%s\t%ld\n", refread.chr, refread.pos);
      while(!endfile)
	{	  
	  if(par.proportion == 1 || (double)rand() / ((double)(RAND_MAX)+(double)(1)) < (double)nreadsleft/(double)(totalreadsleft--))
	    { 
	      nreadsleft--;
	      readsampleline(&read,nextline, &endfile, par.fraglength, par.seqlength, par.qualcutoff, par.paired, par.overlapmode, par.bootstrap); //HAS PUSH IN
	      seqlength = par.paired ? read.pairlength : par.seqlength;
	      fraglength = (par.paired && par.fraglength < read.pairlength) ? read.pairlength : par.fraglength;
	      
	  negpos = set_negpos(&read, &endfile, seqlength, fraglength);
	  //printf("NEG:%ld\t%d\t%ld\t%ld\n", negpos,read.changechr, fragpos, read.pos);
	  //read.changechr = 0;	  
	  set_fragpos(read, &fragpos, &fragend, negpos, fraglength);
	  //printf("FRG:%ld\n", fragpos);
	  reference_reads(read, &refread, &refendfile, par, &nsavedrefchrs, &refnegpos, negpos, &fragpos);
	  while(fragpos <= negpos)
	    {
	      struct countend tempcountend = cntequalhead(fragpos, &pvepos); //CODE	      
	      //printf("TCE1:%f\t%f\t%ld\n", read.count, read.negcount, tempcountend.end);
	      read.count = tempcountend.count + ((fragpos == (read.pos + seqlength - fraglength)) ? read.negcount : 0);
	      if(read.count > par.threshold){ read.count = par.threshold; }
	      if(tempcountend.count){fragend = tempcountend.end; read.prevpos = fragpos;}	      
	      else{ if(fragend < fragpos + fraglength - 1){ fragend = fragpos + fraglength - 1;} }

	      //read.refcount = cntequalhead(fragpos, &refpvepos).count + cntequalhead(fragpos, &refnvepos).count;
	      tempcountend = cntequalhead(fragpos, &refpvepos); read.refcount = tempcountend.count;
	      //if(!empty(&refpvepos)) { printf ("H+: %ld\n", head(&refpvepos).startend[0]);}
	      //if(!empty(&refnvepos)) { printf ("H-: %ld\n", head(&refnvepos).startend[0]);}
	      //printf("R+:%f ", tempcountend.count);
	      
	      //if(tempcountend.count){fragend = tempcountend.end;}
	      tempcountend = cntequalhead(fragpos, &refnvepos); read.refcount += tempcountend.count;
	      //printf("R-:%f\n", tempcountend.count);
	      if(read.refcount > par.refthreshold) { read.refcount = par.refthreshold; }
	      //if(tempcountend.count){fragend = tempcountend.end;}
	      //printf("qar %s %s\n", read.chr, read.prevchr);
	      if(strcmp(read.chr,read.prevchr)) { ref_endchr(read, &curr, par, fragpos, fragend, &lastregion); possend = curr;}  //? 

	      //printf("TCE2:%f\t%f\t%f\t%ld\t%ld\n", read.count, read.negcount, read.refcount, fragpos, fragend);
	      if(read.count > 0 || read.refcount > 0)
		{ 
		  if(strcmp(read.chr,read.prevchr)) 
		    { 
		      //printf("EC %s %s\n", read.chr, read.prevchr);
		      //ref_endchr(read, &curr, par, fragpos, fragend, &lastregion); possend = curr; 
		    }
		  else 
		    {
		      update_total(1, &curr, &possend, negpos, fragpos, fragend, prevfragpos, par, read, &lastregion); 
		      //if(curr.total >= possend.total){ possend = curr;}
		    }
		  while(!empty(&pvepos) && head(&pvepos).startend[0] == fragpos)
		    {  
		      //printf("SUS:%ld\t%ld\t%ld\n", head(&pvepos).startend[0], fragpos,head(&pvepos).startend[1]);
		      sortunshift(fragpos, head(&pvepos).startend[1], head(&pvepos).count, &backpvepos); //CODE

		      push_peakheight(1, fragpos, head(&pvepos).startend[1], head(&pvepos).count, &curr, &possend);

		      shift(&pvepos);
		    }

		  
		  if(par.wig && read.count > 0){push_wigcounts(read, fragpos, prevfragpos, fragend, fraglength, conversion, par.contigrows, wigchr, &wigstart,par.paired, &wigonevalue);}
		  while(!empty(&refpvepos) && head(&refpvepos).startend[0] == fragpos)
		    {
		      sortunshift(fragpos, head(&refpvepos).startend[1], head(&refpvepos).count, &refbackpos);
		      push_peakheight(1, fragpos, head(&refpvepos).startend[1], -1 * par.refpen * head(&refpvepos).count, &curr, &possend);
		      shift(&refpvepos);
		    }
		  while(!empty(&refnvepos) && head(&refnvepos).startend[0] == fragpos)
		    {
		      sortunshift(fragpos, head(&refnvepos).startend[1], head(&refnvepos).count, &refbackpos);
		      push_peakheight(1, fragpos, head(&refnvepos).startend[1], -1 * par.refpen * head(&refnvepos).count, &curr, &possend);
		      shift(&refnvepos);
		    }
		  
		  
		  if(read.count > 0 || strcmp(read.chr,read.prevchr)) { prevfragpos = fragpos; }
		  //printf("y\n");
		}
	      //printf("z\n");
	      strcpy(read.prevchr, read.chr);	
	      
	      //if(!par.paired && read.negcount && (fragpos == (read.pos + seqlength - fraglength))) 
	      if(read.negcount && (fragpos == (read.pos + seqlength - fraglength))) 
		{		
		  /*int i;
		    for(i = 0; i < read.negcount; i++) 
		    {
		      unshift(fragpos, fragpos + fraglength - 1, 1, &backpos); 
		      // printf("UF:%ld\n",fragpos);
		    }
		  */
		  sortunshift(fragpos, fragpos + fraglength - 1, read.negcount, &backnvepos);
		  push_peakheight(1, fragpos, fragpos + fraglength - 1, read.negcount, &curr, &possend);
		}
	      update_fragpos(&fragpos, seqlength, fraglength, read);
	      //printf("FRG2:%ld\n", fragpos);
	      update_fragpos_ref(&fragpos, read.chr, refread.chr);
	      //printf("FRG3:%ld\t%s\t%s\n", fragpos, read.chr, refread.chr);
	    }
	  read.changechr = 0;
	    }
	  else
	    {
	      	if(egets (nextline, 1000, in_fp) == NULL)
		  {endfile = 1;}
	    }
	}
      if(possend.total > 0)
	{ store_results(read, possend, par, &lastregion); }
      
    }
	  
  else
    {
      while(!endfile)
	{
	  if(par.proportion == 1 || (double)rand() / ((double)(RAND_MAX)+(double)(1)) < (double)nreadsleft/(double)(totalreadsleft--))
	    { 	      
	      nreadsleft --;
	      readsampleline(&read, nextline, &endfile, par.fraglength, par.seqlength, par.qualcutoff, par.paired, par.overlapmode, par.bootstrap);
	      seqlength = par.paired ? read.pairlength : par.seqlength;
	      fraglength = (par.paired && par.fraglength < read.pairlength) ? read.pairlength : par.fraglength;
	  //printf("%f\t%f\t%d\t%d\n",(double)rand() / ((double)(RAND_MAX)+(double)(1)), (double)nreadsleft/(double)(totalreadsleft),nreadsleft,totalreadsleft);	 
	      //negpos = par.paired ? set_negpos(&read, &endfile, read.pairlength, (par.fraglength > read.pairlength) ? par.fraglength : read.pairlength) : set_negpos(&read, &endfile, par.seqlength, par.fraglength);
	      negpos = set_negpos(&read, &endfile, seqlength, fraglength);
	  //read.changechr = 0;
	  //printf("NEG:%ld\t%d\t%ld\t%ld\n", negpos,read.changechr, fragpos, read.pos);
	  set_fragpos(read, &fragpos, &fragend, negpos, fraglength); 
	  //set_fragpos(read, &fragpos, &fragend, negpos, (par.paired && par.fraglength < read.pairlength) ? read.pairlength : par.fraglength);
	      //printf("FRG:%ld\n", fragpos);
	  while(fragpos <= negpos)
	    {
	      struct countend tempcountend = cntequalhead(fragpos, &pvepos); //CODE
	      
	      //printf("TCE1:%f\t%f\t%ld\t%f\n", read.count, read.negcount, tempcountend.end, tempcountend.count);
	      read.count = tempcountend.count + ((fragpos == (read.pos + seqlength - fraglength)) ? read.negcount : 0);
 	      if(read.count > par.threshold){ read.count = par.threshold; }
	      if(tempcountend.count > 0){fragend = tempcountend.end; read.prevpos = fragpos;}	      
	      else{ if(fragend < fragpos + fraglength -1){ fragend = fragpos + fraglength - 1;} }
	      //printf("TCE:%f\t%f\t%ld\n", read.count, read.negcount, fragend);
	      if(read.count > 0)
		{
		  //curr.pos++;
		  if(strcmp(read.chr,read.prevchr)) 
		    { 
		      //printf("EC");
		      end_chr(read, &curr, par, fragpos, fragend, &lastregion); possend = curr; 
		      //printf(":%d\n", cnt(&backpos));
		    }
		  else
		    { 
		      update_total(1, &curr, &possend, negpos, fragpos, fragend, prevfragpos, par, read, &lastregion); 
		      //if(curr.pos == 1 || curr.total >= possend.total) { possend = curr; }
		    }		  
		  while(!empty(&pvepos) && head(&pvepos).startend[0] == fragpos)
		    {  
		      //printf("SUS:%s\t%ld\t%f\n", read.chr, head(&pvepos).startend[0], head(&pvepos).count);
		      sortunshift(fragpos, head(&pvepos).startend[1], head(&pvepos).count, &backpvepos); //CODE

		      push_peakheight(1, fragpos, head(&pvepos).startend[1], head(&pvepos).count, &curr, &possend);

		      shift(&pvepos);
		      ///if(par.wig){push_wigcounts(read, fragpos, prevfragpos, fragend, par.fraglength, conversion, par.contigrows, wigchr, &wigstart,par.paired);}
		    }
		  if(par.wig){push_wigcounts(read, fragpos, prevfragpos, fragend, fraglength, conversion, par.contigrows, wigchr, &wigstart,par.paired, &wigonevalue);}

		  strcpy(read.prevchr, read.chr);
		  prevfragpos = fragpos;
		}
	      
	      //if(!par.paired && read.negcount && (fragpos == (read.pos + par.seqlength - par.fraglength))) {
	      if(read.negcount && (fragpos == (read.pos + seqlength - fraglength))) 
		{
		  sortunshift(fragpos, fragpos + fraglength - 1, read.negcount, &backnvepos); 
		  push_peakheight(1, fragpos, fragpos + fraglength - 1, read.negcount, &curr, &possend);
		}
		//int i; for(i = 0; i < read.negcount; i++) {unshift(fragpos, fragpos + par.fraglength - 1, 1, &backpos);}}
	      update_fragpos(&fragpos, seqlength, fraglength, read);
	      //printf("FRG2:%ld\t%ld\n", fragpos,negpos);
	    }
	  read.changechr = 0;
	    }
	  else
	    {
	      	if(egets (nextline, 1000, in_fp) == NULL)
		  {endfile = 1;}
	    }
	}
      if(possend.total>0)
	{ store_results(read, possend, par, &lastregion); }
    } 

  while(!empty(&backpvepos)) 
    { 
      shift(&backpvepos);
      //printf("SH: %s\t%ld\t%f\n", read.chr, head(&backpos).startend[0], head(&backpos).count);
    }
  while(!empty(&backnvepos)) { shift(&backnvepos); }
  while(!empty(&refbackpos)) { shift(&refbackpos); }

  if(!empty(&wigcounts)) {export_wiggle(conversion, par.contigrows, wigchr, &wigstart, par.fraglength, par.paired, wigonevalue);}
  if(par.comp) { 
    compare_lastregion(&lastregion); 
    print_out_overlap();
  }
  else 
    {
      if((lastregion.count - lastregion.refcount*par.refpen) >= par.min_above_bg && lastregion.total >= par.resultcutoff)
	{ print_out_lastregion(&lastregion, par.binom_p, par.newpeakversion); }
    }

  while(par.ref && !refendfile)
    {
      readrefline(&refread, &refendfile, par.reffraglength, par.paired, par.bootstrap);
      if(strcmp(refread.chr, refread.prevchr))
	{ printf("Discarding ref. chromosome %s !\n", refread.chr);}
      //strcpy(refread.prevchr, refread.chr);
    }
  return(0);
}



void store_results(struct readinfo read, struct region possend, struct param par, struct region *lastregion)
 {
   struct region pvestrandregion, olr;
   long int fragpos, fragstart, prevfragpos;
   //printf("SR:%d\t%f\t%f\t%f\n", cnt(&backpvepos) + cnt(&backnvepos), (possend.count - possend.refcount), possend.total, possend.peakheightmax);
   //if(((possend.count - possend.refcount*par.refpen) >= par.min_above_bg && possend.total >= par.resultcutoff) || (strcmp(read.chr,read.prevchr)))
   if(((possend.count - possend.refcount*par.refpen) >= par.min_above_bg && possend.total >= par.resultcutoff) || (strcmp(read.chr,read.prevchr)))
     {
       strcpy(read.chr, read.prevchr);
       //if(num_above_bg >= par.min_above_bg && possend.total >= par.resultcutoff)
       // { pvestrandregion = possend; }

       if(!empty(&peakheights)) 
	 { 
	   calculate_peak_summit(&possend, head(&peakheights).startend[0], tail(&peakheights).startend[0], &possend);   
	 }


       /*
	  int i;
	  for(i = 0; i < cnt(&peaksummits); i++)
	   {
	 printf("%ld ", (long int)return_position(i, &peaksummits));
	   }
	  printf("\n");
       */

       if(!empty(&peaksummits))
	 {
	   possend.summit = median(&peaksummits);
	 }

       

       pvestrandregion = possend;
       if((possend.count - possend.refcount*par.refpen) < par.min_above_bg || possend.total < par.resultcutoff)
	 {pvestrandregion.total = 0; pvestrandregion.print_total = 0;}
       else
	 {
	   //fprintf(out_fp,"%s\t%ld\t%ld\t%d\t%ld\t%d\t%f\t%d\t+\n", possend.chr, possend.start, possend.end, possend.count, possend.end-possend.start+1, possend.pos, possend.total,possend.refcount);
	   // fprintf(out_fp,"%s\t%ld\t%ld\t%f\t%ld\t%d\t%f\t%f\t%f\t%f\t%ld\t%ld\t+\n", possend.chr, possend.start, possend.end, possend.count, possend.end-possend.start+1, possend.pos, possend.total, possend.refcount, possend.peakheightmax, possend.summit, possend.seen_start, possend.seen_end);
	   // printf("%s\t%ld\t%ld\t%f\t%ld\t%d\t%f\t%f\t%f\t%f\t%ld\t%ld\t+\n", possend.chr, possend.start, possend.end, possend.count, possend.end-possend.start+1, possend.pos, possend.total, possend.refcount, possend.peakheightmax, possend.summit, possend.seen_start, possend.seen_end);
	   
	   // while(!empty(&peaksummits)) { printf("%f ", shift(&peaksummits).count); } printf("\n");
	 }

       int backstrand = set_back_fragpos(&fragpos, &fragstart);

       if(!par.paired && backstrand == 1)
	 { read.pos = fragpos - par.fraglength + par.seqlength ; }
       else{ read.pos = fragpos; }
       
       olr = pvestrandregion;
       //save fragpos :remove
       read.count = 0;
       read.refcount = 0;
       
       
       
       //printf("BFRG:%ld\t%ld\n", fragpos, fragstart);
       struct region back = new_region(read, par, fragpos, fragstart, -1);
       //strcpy(back.chr, possend.chr);
       possend = back;

       prevfragpos = fragpos;
       while(!empty(&backpvepos) || !empty(&backnvepos))
	 {
	   backstrand = set_back_fragpos(&fragpos, &fragstart);
	   //printf("BFRG2:%ld\t%ld\n", fragpos, fragstart);
	   if(!par.paired && backstrand == 1)
	     { read.pos = fragpos - par.fraglength + par.seqlength ; }
	   else{ read.pos = fragpos; }
	   
	   
	   if(!empty(&refbackpos) && head(&refbackpos).startend[1] > fragpos) 
	     { fragpos = head(&refbackpos).startend[1]; }

	   //read.count = sample_back_count(fragpos, &fragstart, par.threshold, &backpos);
	   //read.count = sample_back_count(fragpos, &fragstart, par.threshold, &backpvepos) + 
	   //           sample_back_count(fragpos, &fragstart, par.threshold, &backnvepos);
	   read.negcount = sample_back_count(fragpos, &fragstart, par.threshold, &backpvepos);
	   read.count = read.negcount + sample_back_count(fragpos, &fragstart, par.threshold, &backnvepos);
	   read.refcount = ref_back_count(fragpos, par.refthreshold, &refbackpos);
 
	   //printf("SBC:%ld\t%ld\t%f\t%f\n", fragpos, fragstart, read.count, read.refcount);
	   //if(read.count > 0){ back.pos++; }
	   //update_total(-1, &back, &possend, fragpos+1-backstrand, fragpos, fragstart, prevfragpos, par, read, lastregion); //negpos?? 
	   update_total(-1, &back, &possend, fragpos, fragpos, fragstart, prevfragpos, par, read, lastregion); //negpos?? 

	   //if(back.pos == 1 || back.total >= possend.total) { possend = back; }
	   if(back.total <= 0.00001)
	     { 
	       output_regions(possend, &olr, pvestrandregion, par);
	       back = new_region(read, par, fragpos, fragstart, -1);
	       possend = back;
	     }	   

	   while(!empty(&backpvepos) && head(&backpvepos).startend[1] == fragpos)
	     {
	       push_peakheight(-1, fragpos, head(&backpvepos).startend[0], head(&backpvepos).count, &back, &possend);
	       //printf("SH+: %s\t%ld\t%f\n", read.chr, head(&backpvepos).startend[0], head(&backpvepos).count);
	       shift(&backpvepos);
	     }

	   while(!empty(&backnvepos) && head(&backnvepos).startend[1] == fragpos)
	     {
	       push_peakheight(-1, fragpos, head(&backnvepos).startend[0], head(&backnvepos).count, &back, &possend);
	       //printf("SH-: %s\t%ld\t%f\n", read.chr, head(&backnvepos).startend[0], head(&backnvepos).count);
	       shift(&backnvepos);
	     }

	   
	   while(!empty(&refbackpos) && head(&refbackpos).startend[1] == fragpos)
	     {
	       push_peakheight(-1, fragpos, head(&refbackpos).startend[0], -1 * par.refpen * head(&refbackpos).count, &back, &possend);
	       shift(&refbackpos);
	     }
	   //CODE store_results || output_regionssss
	   if(read.count > 0) { prevfragpos = fragpos; }
	 }
       
       //printf("ORSLAST:%f\t%f\t%f\n", possend.count, possend.refcount, possend.total);
       if((possend.count - possend.refcount*par.refpen) >= par.min_above_bg && possend.total >= par.resultcutoff)
	 { output_regions(possend, &olr, pvestrandregion, par); }

       //printf("CPR\n");
       *lastregion = check_and_print_regions(olr, &(*lastregion), par.comp, par.binom_p, par.newpeakversion);
       
       //printf("%s\t%ld\t%ld\t%ld\t%d\t%f\t%f\t%ld\n", (*lastregion).chr, (*lastregion).start, (*lastregion).end, (*lastregion).end-(*lastregion).start+1, (*lastregion).pos, (*lastregion).print_total, (*lastregion).peakheightmax, (long int)floor((*lastregion).summit));

       while(!empty(&refbackpos))
	 { shift(&refbackpos); }
     }
 }
	     
