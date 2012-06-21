#include "SWEMBL.h"

void reference_reads(struct readinfo read, struct readinfo *refread, int *refendfile, struct param par, int *nsavedrefchrs, long int *refnegpos, long int negpos, long int *fragpos)
  {
    //printf("%s\t%s\t%s\t%s\n", read.chr, read.prevchr, (*refread).chr, (*refread).prevchr);
    unsigned int refseqlength=0;
    int reffraglength=0;
    
    if(strcmp(read.chr, read.prevchr))
      {
	//printf("Del\n");
	while(!empty(&refpvepos))
	  { shift(&refpvepos); }
	while(!empty(&refnvepos))
	  { shift(&refnvepos); }
      }
    
    //printf("B\n");
    while(strcmp(read.chr, read.prevchr) && !strcmp((*refread).chr, read.prevchr) && !(*refendfile)) //End of sample reads on chr, ignore later ref reads
      {
	//printf("G\n");
	/*int ms;
	  int exists_savedrefs_chr = 0;
	  for(ms = 0; ms < maxsavedrefchrs; ms++)
	  {
	  printf("H\n");
	  if(!strcmp(savedrefs[ms].chr, chr))
	  { exists_savedrefs_chr = 1; break;}
	  }
	  if(exists_savedrefs_chr)
	  {
	  while(!empty(&savedrefpos))
	  { shift(&savedrefpos);}
	  strcpy(refchr, chr);
	  }
	*/
	
	if(!empty(&savedrefpos))
	  {
	    //printf("AC\t%s\t%s\t%s\n",read.chr,read.prevchr,(*refread).chr);
	    while(!empty(&savedrefpos))
	      { shift(&savedrefpos); }
	    strcpy((*refread).chr, read.chr);
	  }
	else
	  {
	    readrefline(refread, refendfile, par.reffraglength, par.paired, par.bootstrap); 
	    refseqlength = par.paired ? (*refread).pairlength : par.refseqlength;
	    reffraglength = (par.paired && par.reffraglength < (*refread).pairlength) ? (*refread).pairlength : par.reffraglength;
	    
	    //printf("AB\t%s\t%s\t%s\t%s\n",read.chr,read.prevchr,(*refread).chr, (*refread).prevchr);
	    //*refnegpos = (*refread).pos - par.reffraglength + par.refseqlength;  
	    *refnegpos = 0;
	  }
      }
    

    if(strcmp(read.chr,read.prevchr))
      {
	//printf("YY\n");
	int ms = 0;
	int exists_savedrefs_chr = 0;
	//for(ms = 0; ms < par.maxsavedrefchrs; ms++)
	//for(ms = 0; ms < *nsavedrefchrs+1; ms++)
	//  {
	//  if(savedrefs[ms].chr && !strcmp(savedrefs[ms].chr, read.chr))
	//    { exists_savedrefs_chr = ms+1; break;}
	  
	//}


	if(!emptychr(&refchrs_pve)) 
	  {	    	
	    chrelem *c_pve, *c_nve;
	    c_pve = (&refchrs_pve) -> head;
	    c_nve = (&refchrs_nve) -> head;
	    if(strcmp(c_pve->chr, read.chr)) 
	      {
		for(ms = 1; ms < chrcnt(&refchrs_pve) ; ms++)
		  {
		    //printf("ms %d\t%s\t%s\n", ms, c_pve->next->chr, read.chr);
		    if(!strcmp(c_pve->next->chr, read.chr))
		      {exists_savedrefs_chr = ms + 1; break;}
		    c_pve = c_pve->next;
		    c_nve = c_nve->next;
		  }
	      }
	    else
	      { exists_savedrefs_chr = 1; }
	  
	  


	      if(exists_savedrefs_chr) //New chromosome, reference reads previously saved
		{
		  while(!empty(&savedrefpos))
		    { shift(&savedrefpos);}
		  //printf("Ea %d\n", chrcnt(&refchrs_pve));
		  //refpvepos = (*c -> next -> stk);
		  if(exists_savedrefs_chr == 1)
		    {
		      //savedpverefs = (*c_pve -> stk);
		      assign(&savedpverefs, c_pve -> stk);
		      ((&refchrs_pve) -> cnt)--;
		      if(!emptychr(&refchrs_pve)) { (&refchrs_pve) -> head = (&refchrs_pve) -> head -> next; }
		      //shiftchr(&refchrs_pve);
		      //savednegrefs = (*c_nve -> stk);
		      assign(&savednegrefs, c_nve -> stk);
		      //shiftchr(&refchrs_nve);
		      ((&refchrs_nve) -> cnt)--;
		      if(!emptychr(&refchrs_nve)) { (&refchrs_nve) -> head = (&refchrs_nve) -> head -> next; }
		    }
		  else
		    {
		      if(exists_savedrefs_chr == chrcnt(&refchrs_pve))
			{
			  //savedpverefs = c -> next -> stk;
			  assign(&savedpverefs, c_pve->next->stk);
			  (&refchrs_pve) -> tail = c_pve;
			  ((&refchrs_pve) -> cnt)--;
		      //printf("Tla %d %s\n", cnt(c_pve->next->stk), c_pve->next->chr);
			  assign(&savednegrefs, c_nve->next->stk);
			  (&refchrs_nve) -> tail = c_nve;
			  ((&refchrs_nve) -> cnt)--;
			  //printf("Tlb %d %s\n", cnt(&savednegrefs), c_nve->next->chr);
			}
		      else
			{
			  //savedpverefs = (*c_pve -> next -> stk);
			  assign(&savedpverefs, c_pve->next->stk);
			  c_pve -> next = c_pve -> next -> next;
			  ((&refchrs_pve) -> cnt)--;
			  //savednegrefs = (*c_nve -> next -> stk);
			  assign(&savednegrefs, c_nve->next->stk);
			  c_nve -> next = c_nve -> next -> next;
			  ((&refchrs_nve) -> cnt)--;
			}
		    }
		  (*nsavedrefchrs)--;
	    

		  while(cnt(&savedpverefs) || cnt(&savednegrefs))
		    {
		      if(empty(&savednegrefs))
			{ 
			  push(head(&savedpverefs).startend[0],head(&savedpverefs).startend[1], head(&savedpverefs).count, &savedrefpos); shift(&savedpverefs);
			}
		      else
			{ 
			  if(empty(&savedpverefs) || head(&savednegrefs).startend[0] < head(&savedpverefs).startend[0])
			    { 
			      push(head(&savednegrefs).startend[0],head(&savednegrefs).startend[1], head(&savednegrefs).count, &savedrefpos); shift(&savednegrefs); 
			    } 
			  else 
			    { 
			      push(head(&savedpverefs).startend[0],head(&savedpverefs).startend[1], head(&savedpverefs).count, &savedrefpos); shift(&savedpverefs);
			    } 
			}
		    }
		 
		  strcpy((*refread).chr, read.chr);
		  strcpy((*refread).prevchr, (*refread).chr); 
		  //printf("C:%d\n",cnt(&savedrefpos));
		  ///(*refread).pos = shift(&savedrefpos).startend[0];
		    (*refread).pos = head(&savedrefpos).startend[0];
		    //printf("D:%s\n",(*refread).chr);
		    *refnegpos = (*refread).pos - reffraglength + refseqlength;
		    strcpy((*refread).strand, pvestring);

		}
	  }
	if(!par.quietmode){printf("%s\tREF:%s\n", read.chr, (*refread).chr);}
      
      }
    
    
	///int rs;
	  ///int exists_savedrefs_refchr = 0;
	  ///for(rs = 0; rs < *nsavedrefchrs+1; rs++)
	  ///{
	  ///if(!strcmp(savedrefs[rs].chr, (*refread).chr))
	  ///{ exists_savedrefs_refchr = rs+1; break;}
	  ///
	  ///

    int ms = 0;
    int exists_savedrefs_refchr = 0;
    chrelem *c_pve, *c_nve;
    //c = malloc(sizeof(chrelem);
    if(!emptychr(&refchrs_pve)) 
      {
	c_pve = (&refchrs_pve) -> head;
	c_nve = (&refchrs_nve) -> head;
	//c_pve = firstchr(&refchrs_pve);
	
	if(strcmp(c_pve->chr, read.chr)) 
	  {
	    for(ms = 1; ms < chrcnt(&refchrs_pve) - 1; ms++)
	      {
		if(!strcmp(c_pve->next->chr, read.chr))
		  {exists_savedrefs_refchr = ms + 1; break;}
		c_pve = c_pve->next;
		c_nve = c_nve->next;
	      }
	  }
	else
	  { exists_savedrefs_refchr = 1; }
      }

    //printf("EX%d\t%d\n",exists_savedrefs_refchr,*nsavedrefchrs);
    
    ///if((exists_savedrefs_refchr) || *nsavedrefchrs < par.maxsavedrefchrs)
	  ///  {
    while(strcmp(read.chr,read.prevchr) && strcmp((*refread).chr,read.chr) && !(*refendfile)) //New chr, different ref chr
      {
	if(strcmp((*refread).chr, (*refread).prevchr))
	  {
		///if(*nsavedrefchrs == par.maxsavedrefchrs)
		///  {  //readrefline(refread, refendfile, par.reffraglength, par.paired); 
		    //printf("Break\n");
		///    break;}
		///else
		///  { 
		    ///strcpy(savedrefs[++(*nsavedrefchrs) -1].chr, (*refread).chr);
		    ///savedrefs[*nsavedrefchrs -1].pvecount = 0;
		    ///savedrefs[*nsavedrefchrs -1].negcount = 0;
	    addchr((*refread).chr, &refchrs_pve);
	    addchr((*refread).chr, &refchrs_nve);
	    (*nsavedrefchrs)++;
	    
		    ///printf("Saving: %s\t%d\n", savedrefs[*nsavedrefchrs - 1].chr, *nsavedrefchrs);
	    printf("Saving: %s\t%d\n", (*refread).chr, *nsavedrefchrs);
		      /// }
	  }
	if((*refread).qual >= par.refqualcutoff)
	  {
		//if(strcmp((*refread).strand, negstring) && (!par.paired || (*refread).pairflag == 18 ) )
	    
	    if((*refread).pvecount != 0)
	      {
		 ///push((*refread).pos, (*refread).pos + (*refread).pairlength - 1, (*refread).pvecount, &savedpverefs); 
		     ///savedrefs[*nsavedrefchrs - 1].pvecount++;
		    //push((*refread).pos-par.reffraglength+par.refseqlength, &savedpverefs); 
		pushchr((*refread).pos, (*refread).pos + (*refread).pairlength - 1, (*refread).pvecount, &refchrs_pve);
	      }
		//else
		//{
		  //   if(strcmp((*refread).strand, pvestring) && !par.paired)
		//     {
	     if((*refread).negcount)
	       {
		 ///push((*refread).pos - par.reffraglength + par.refseqlength, (*refread).pos+par.refseqlength-1, (*refread).negcount, &savednegrefs);
		       ///savedrefs[*nsavedrefchrs - 1].negcount++;
		       //push((*refread).pos, &savednegrefs);
		       //}
		 pushchr((*refread).pos - reffraglength + refseqlength, (*refread).pos+refseqlength-1, (*refread).negcount, &refchrs_nve);
	       }
	   
	    

	  }
	    //printf("AM:%s\t%s\n",(*refread).chr,(*refread).prevchr);
	    
	readrefline(refread, refendfile, par.reffraglength, par.paired, par.bootstrap);
	refseqlength = par.paired ? (*refread).pairlength : par.refseqlength;
	reffraglength = (par.paired && par.reffraglength < (*refread).pairlength) ? (*refread).pairlength : par.reffraglength;
      }
	// printf("Saved: %s\t%d\t%ld\t%ld\n", savedrefs[*nsavedrefchrs - 1].chr, *nsavedrefchrs,savedrefs[*nsavedrefchrs - 1].pvecount,savedrefs[*nsavedrefchrs - 1].negcount);
	///}
    
    ///int ms;
      ///int exists_savedrefs_chr = 0;
	  /// for(ms = 0; ms < *nsavedrefchrs+1; ms++)
      /// {
	  ///if(!strcmp(savedrefs[ms].chr, read.chr))
	   ///  { exists_savedrefs_chr = ms+1; break;}
      ///}

    int exists_savedrefs_chr = 0;
    if(!emptychr(&refchrs_pve)) 
      {
	c_pve = (&refchrs_pve) -> head;
	if(strcmp(c_pve->chr, read.chr)) 
	  {
	    for(ms = 1; ms < chrcnt(&refchrs_pve) - 1; ms++)
	      {
		if(!strcmp(c_pve->next->chr, read.chr))
		  {exists_savedrefs_chr = ms + 1; break;}
		c_pve = c_pve->next;
	      }
	  }
	else
	  { exists_savedrefs_chr = 1; }     
      }

    //free(c);
    
      //printf("RNC:%ld\t%s\t%d\n", *refnegpos, (*refread).chr,exists_savedrefs_chr);
    while((*refnegpos <= negpos) && !strcmp((*refread).chr,read.chr) && (!(*refendfile) || !empty(&savedrefpos))) //Within a chromosome
      {
	//printf("D:%s\t%ld\t%ld\n",(*refread).chr,*refnegpos,negpos);
	//if(exists_savedrefs_chr)
	//{

	//printf("E\n");
	//push((*refread).pos,(*refread).pos + (*refread).pairlength -1, &refpvepos);
	if(!empty(&savedrefpos) && (*refread).count)
	  {	    
	    push((*refread).pos,(*refread).pos + (*refread).pairlength -1, (*refread).count, &refpvepos);
	    (*refread).pos = shift(&savedrefpos).startend[0];
	    *refnegpos = (*refread).pos - reffraglength + refseqlength;
	    //printf("S:%ld\n",*refnegpos);
	  }
	//else { strcpy(refchr,"0"); }
		   
	//}
	else
	  {
	    //printf("A\t%d\t%s\n", (*refread).qual, (*refread).strand);
	    if((*refread).qual >= par.refqualcutoff)
	      {
		//if(strcmp((*refread).strand, negstring) && (!par.paired || (*refread).pairflag == 18 ) )
		if((*refread).pvecount)
		  { push((*refread).pos, (*refread).pos + (*refread).pairlength-1, (*refread).pvecount, &refpvepos); }
		//else
		// { 
		//    if(strcmp((*refread).strand, pvestring) && !par.paired)
		if((*refread).negcount)
		      { push((*refread).pos - reffraglength + refseqlength, (*refread).pos+refseqlength-1, (*refread).negcount, &refnvepos); }
		      //  }
	      }
	    readrefline(refread, refendfile, par.reffraglength, par.paired, par.bootstrap);
	    refseqlength = par.paired ? (*refread).pairlength : par.refseqlength;
	    reffraglength = (par.paired && par.reffraglength < (*refread).pairlength) ? (*refread).pairlength : par.reffraglength;
	    *refnegpos = (*refread).pos - reffraglength + refseqlength;
	  }
      }
   
    (*refread).changechr = 0;     ///???
    
    
    if(!empty(&pvepos) && head(&pvepos).startend[0] < negpos)
      { *fragpos = head(&pvepos).startend[0]; }
    else
      { *fragpos = negpos; }
    
    if(!empty(&refpvepos) && head(&refpvepos).startend[0] < *fragpos)
      { *fragpos = head(&refpvepos).startend[0]; }
    if(!empty(&refnvepos) && head(&refnvepos).startend[0] < *fragpos)
      { *fragpos = head(&refnvepos).startend[0]; }
    
      
    //printf("Negpos\t%ld\n",negpos);
    /*
      if(!empty(&pvepos)) 
      { printf("Pvepos\t%ld\n",head(&pvepos)); }
      if(!empty(&refnvepos)) 
      { printf("Refnvepos\t%ld\n",head(&refnvepos)); } 
      if(!empty(&refpvepos)) 
      { printf("Refpvepos\t%ld\n",head(&refpvepos)); }
      */
    
    //printf("A2%ld\t%ld\n",*fragpos,negpos);
  }


/*int refcts(long int fragpos, int threshold)
  {
    int refcount = 0;
    while(!empty(&refpvepos) && head(&refpvepos) == fragpos)
      {
	refcount++;
	shift(&refpvepos);
      }
    
    while(!empty(&refnvepos) && head(&refnvepos) == fragpos)
      {
	refcount++;
	shift(&refnvepos);
      }
    
    if(refcount >= refthreshold)
      { refcount = refthreshold; }
    / * * /printf("%ld\t%d\t%d\t%d\t%d\t%s\t%s\n", fragpos, pvecount, negcount, count, refcount, chr, prevchr);
    return(refcount);
3  }
*/


void ref_endchr(struct readinfo read, struct region *curr, struct param par, int fragpos, int fragend, struct region *lastregion)
   {
     //printf("%s\n", read.chr);
     //printf("REC:%s\t%s\n", read.chr, read.prevchr);
     if(strcmp(read.prevchr,""))
       {
	 //genomelength++;
	 //pvestrandregion.hightot=0;
	 //strcpy(read.chr, read.prevchr);
	 store_results(read, *curr, par, lastregion); 
	while(!empty(&backpvepos))
	  { shift(&backpvepos); }
	while(!empty(&backnvepos))
	  { shift(&backnvepos); }
	 while(!empty(&refbackpos))
	 { shift(&refbackpos); }
 
       }
     
     //else { strcpy(read.prevchr, read.chr) ; }  // ???
     *curr = new_region(read, par, fragpos, fragend, 1);
   }




void update_fragpos_ref(long int *fragpos, char readchr[], char refreadchr[])
{ 
  if(!empty(&refpvepos) && !strcmp(readchr,refreadchr) && head(&refpvepos).startend[0] < *fragpos)
    { *fragpos = head(&refpvepos).startend[0]; }
  if(!empty(&refnvepos) && !strcmp(readchr,refreadchr) && head(&refnvepos).startend[0] < *fragpos)
    { *fragpos = head(&refnvepos).startend[0]; }
}


double ref_back_count(long int fragpos, double refthreshold, stack *stk)
{
  double refcount = 0;
  /*  
   while(!empty(&refbackpos) && head(&refbackpos).startend[1] == fragpos)
     {
       refcount += head(&refbackpos).count;
       shift(&refbackpos);
     }
  */
  
  if(!empty(stk))
    {
      elem *p = stk -> head;
      int i = 0;
     
      if(p -> d[1] == fragpos) 
	{ 
	  refcount += p -> s; 
	  while(i < (cnt(stk) -1)  && p -> next->d[1] == fragpos)
	    {
	      p = p->next;	      
	      refcount += p -> s;
	      i++;
	    }
	}      
      // printf("SBC\t%d\t%f\n", i, count);
    }  
  
  if(refcount > refthreshold){ refcount = refthreshold;} //REF THRESHOLD??
  return(refcount);
}


