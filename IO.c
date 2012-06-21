#include "SWEMBL.h"

struct param get_parameters(char infile[1000], char outfile[1000], char reffile[1000], char wigfile[1000], char contigfile[1000], char compfile[1000], char compoverlapfile[1000], char regoverlapfile[1000], int argc, char **argv)
{
  struct param par;
  par.bg = 0;
  par.posbg = 0;
  par.longbg = 0;
  par.min_above_bg = 0;
  par.threshold = 5;
  par.seqlength = 0;
  par.qualcutoff = 0;
  par.resultcutoff = 0;
  par.pen_inc = 70;
  par.maxtotal = -1;
  par.fraglength = 0;
  par.ref = 0;
  par.reffraglength = 0;
  par.refseqlength = 0;
  par.refpen = 1;
  par.refthreshold = 0;
  par.refqualcutoff = 0;
  par.killifnoref = 0;
  par.paired = 0;
  par.wig = 0;
  par.contigrows = 200;
  par.comp = 0;
  overlap.print = 0;
  par.proportion = 1;
  par.totalreads = 0;
  par.totalrefreads = 0;
  par.quietmode = 0;
  par.overlapmode = 0;
  par.bootstrap = 0;
  par.binom_p = 1;
  par.newpeakversion=0;
  /*outfile = '';*/
  
  //unsigned int totalreads = 0;
  double relativebg = 0;

  opterr = 0;
  int c;
  int index;
 
  while((c = getopt (argc, argv, "i:o:r:w:a:b:p:P:m:t:j:l:f:q:c:d:s:x:G:n:N:K:R:u:D:AOCBEMSFhezvygTV")) != -1)
   switch(c)
     {
     case 'i': strcpy(infile,optarg); break;
     case 'o': strcpy(outfile,optarg); break;
     case 'r': strcpy(reffile,optarg); break;
     case 'w': strcpy(wigfile,optarg); break;
     case 'a': strcpy(compfile,optarg); break;
     case 'A': strcpy(filemode,"a"); break;
     case 'b': par.bg = atof(optarg); break;
     case 'p': par.posbg = atof(optarg); break;
     case 'P': par.longbg = atof(optarg); break;
     case 'm': par.min_above_bg = atoi(optarg); break;
     case 't': par.threshold = atof(optarg); break;
     case 'j': par.refthreshold = atof(optarg); break;
     case 'l': par.seqlength = atoi(optarg); break;
     case 'f': par.fraglength = atoi(optarg); break;
     case 'N': par.totalreads = atoi(optarg); break;
     case 'K': par.totalrefreads = atoi(optarg); break;
     case 'R': relativebg = atof(optarg); break;
       /*case 'C': countfile = 1; break;
	 case 'B': bedfile = 1; break;
	 case 'E': elandfile = 1; break;
	 case 'M': maqfile = 1; break;*/
     case 'C': strcpy(filetype,"C"); break;
     case 'B': strcpy(filetype,"B"); break;
     case 'E': strcpy(filetype,"E"); break;
     case 'M': strcpy(filetype,"M"); break;
     case 'S': strcpy(filetype,"S"); break;
     case 'F': strcpy(filetype,"F"); break;
     case 'q': par.qualcutoff = atof(optarg);  break;
     case 'c': par.resultcutoff = atof(optarg); break;
     case 'd': par.pen_inc = atoi(optarg); break;
     case 'n': par.maxtotal = atof(optarg); break;
     case 'x': par.refpen = atof(optarg); break;       
     case 'z': zip = 1; break;
     case 'e': par.paired = 1; break;
     case 'V': par.newpeakversion = 1; break;
     case 's': strcpy(contigfile,optarg); break;
     case 'G': par.contigrows = atoi(optarg); break;
     case 'O': overlap.print = 1; break; 
     case 'u': par.proportion = atof(optarg); break;
     case 'y': par.quietmode = 1; break;
     case 'g': par.overlapmode = 1; break;
     case 'T': par.bootstrap = 1; break;
     case 'D': rseed = atoi(optarg); break;
     case 'v': printf("SWEMBL version 3.6\n"); exit (0);
     case 'h': print_help();
     case '?': 
       if (optopt == 'i')
	 fprintf(stderr, "Option -%c requires an argument.\n", optopt);
       else if (isprint(optopt))
	 fprintf(stderr, "Unknown option '-%c'.\n", optopt);
       else
	 fprintf(stderr,"Unknown option character '\\x%x'.\n", optopt);
       exit(1);
     default: 
       abort();
     }

 for (index = optind; index < argc; index++)
   printf ("Non-option argument %s\n", argv[index]);

 //If Input file is not specified, print help, then exit.

 if(!strcmp(infile,""))
   {
     print_help();
   }

 if(!strcmp(outfile,""))
    {
      strcpy(outfile,infile);
      strcat(outfile, ".SWEMBL.3.6.txt");
    }

 if(strcmp(reffile,""))
   {
     par.ref = 1;
     if(!par.refthreshold) { par.refthreshold = par.threshold; }
   }

 if(strcmp(wigfile,""))
   {
     par.wig = 1;
   }

 if(strcmp(compfile,""))
   {
     par.comp = 1;
   }

 if(!strcmp(filetype,""))
   {
     par.overlapmode ? strcpy(filetype,"B") : strcpy(filetype,"F");
   }

 if(zip && !strcmp(filetype,"F"))
   {
     printf("Cannot open a gzipped BAM file!\n");
     exit(1);
   }

 if(!par.comp)
   {
     overlap.print = 0;
   }

 if(overlap.print)
   {    
     strcpy(compoverlapfile, outfile);
     strcat(compoverlapfile, ".overlapscomp");     
     strcpy(regoverlapfile, outfile);
     strcat(regoverlapfile, ".overlapsregion");
   } 
 
 if(par.proportion > 1 || par.proportion <= 0)
   {
     printf("Invalid proportion entered : %f !\n", par.proportion);
     exit(1);
   }

 par.refpen *= par.proportion;

 if(relativebg && !par.totalreads)
   {
     par.totalreads = count_sample_lines(&par, infile, 1);
   }

 if(par.ref && relativebg && !par.totalrefreads)
   {
     par.totalrefreads = count_sample_lines(&par, reffile, 0);
   }

 if(relativebg && (par.totalreads <= 0 || (par.ref && par.totalrefreads <= 0)))
   {
     printf("Number of reads must be greater than zero.\n");
     exit(1);
   }
 
 if(par.totalreads && relativebg)
   {
     par.posbg = (par.posbg == 0) ? relativebg * par.totalreads * par.proportion / 3000000 : par.posbg;
     par.longbg = (par.longbg == 0) ? relativebg * par.totalreads * par.proportion / 600000 : par.longbg;
     par.maxtotal = (par.maxtotal == -1) ? relativebg * 300 * par.proportion : par.maxtotal;
     par.min_above_bg = (par.min_above_bg == 0) ? (((relativebg * par.totalreads * par.proportion / 7500) > 6) ? ceil(relativebg * par.totalreads *par.proportion / 7500) : 6) : par.min_above_bg;
     if(par.ref)
       {
	 par.refpen *= ((double)par.totalreads /(double)par.totalrefreads);
	 par.refthreshold /= par.refpen;
       }
   }
 else
   {
     par.posbg = (par.posbg == 0) ? 0.01 : par.posbg;
     par.longbg = (par.longbg == 0) ? 0.02 : par.longbg;
     par.maxtotal = (par.maxtotal == -1) ? 0 : par.maxtotal;
     par.min_above_bg = (par.min_above_bg == 0) ? 8 : par.min_above_bg;
   }

 if(par.overlapmode)
   {
     par.posbg = 0.9;
     par.longbg = 100;
     par.min_above_bg = 1; 
     par.paired = 1;
     overlap.mode = 1;
   }

 return(par);
}

void print_help()
{
  printf("\nProgram to call peaks from genomic and other data.\n\nMandatory command line options:\n-i Path to input file\n\nOptional command line options taking no arguments:\n-e Paired end reads\n-V New peak version, extending peaks only as far as evidence permits\n-M MAQ format [default]\n-B BED format\n-E ELAND format\n-S SAM format\n-F BAM format\n-C count format\n-z File is gzipped\n-O Create files of nearest features with comparison file (-a)\n-g Overlap mode for two bed files\n-T Bootstrap input\n-y Quiet mode, only print warnings\n-A Append output file\n-v Print version number and exit\n-h This help page\n\nOptional command line options taking arguments [defaults in square brackets]:\n-o Path to output file [Input file.'.SWEMBL.3.6.txt']\n-r Path to reference (Input) file\n-w Path to wiggle track file\n-a Path to comparison file\n-s Path to coordinate file\n-G Number of rows in coordinate file [200]\n-b Penalty applied to read count [0]\n-p Penalty applied to a gap of one base pair [0.01]\n-l Read length [calculated from first read in file for BED and MAQ formats]\n-f Fragment length for extending fragments [Read length]\n-d Gap at which penalty per base pair increases [70]\n-P Penalty per base pair applied after gap specifed above [0.02]\n-t Threshold value for number of sample reads starting at one base pair [5]\n-j Threshold value for number of reference reads starting at one base pair [0 = -t value]\n-q Mapping quality filter value (MAQ only) [0]\n-m Minimum read count in a peak [8]\n-c Minimum score for a peak [0]\n-n Maximum score threshold [0 = infinite]\n-N Number of input reads [0]\n-R Relative background (depends on -N) [0]\n-K Number of reference reads [0]\n-x Penalty factor for reference reads [1]\n-u Proportion of reads sampled [1]\n-D Random seed (only used for sampling and bootstrapping) [time]\n\n"); exit(1);
}


void open_files(int ref, int wig, int comp, char *infile, char *outfile, char *reffile, char *wigfile, char *contigfile, char *compfile, char *compoverlapfile, char *regoverlapfile)
{
  if(zip)
    {
      if((in_fp = gzopen(infile, "r")) == NULL)
	{
	  printf("\n\n*** %s does not exist ***\n", infile);
	  exit(1);
	}

      if(ref && (ref_fp = gzopen(reffile, "r")) == NULL)
	{
	  printf("\n\n*** %s does not exist ***\n", reffile);
	  exit(1);
	}
    }
  else
    {
      if(strcmp(filetype,"F"))
	{
	  if((in_fp = fopen(infile, "r")) == NULL)
	    {
	      printf("\n\n*** %s does not exist ***\n", infile);
	      exit(1);
	    }
	  
	  if(ref && (ref_fp = fopen(reffile, "r")) == NULL)
	    {
	      printf("\n\n*** %s does not exist ***\n", reffile);
	      exit(1);
	    }
	}

      else
	{
	  strcpy(filetype, "S");
	  char samtoolscommand[1000] = "samtools view ";
	  strcat(samtoolscommand, infile);
      
	  if((in_fp = popen(samtoolscommand, "r")) == NULL)
	    {
	      printf("\n\n*** %s does not exist ***\n", infile);
	      exit(1);
	    }

	  strcpy(samtoolscommand, "samtools view ");
	  strcat(samtoolscommand, reffile);
	  
	  if(ref && (ref_fp = popen(samtoolscommand, "r")) == NULL)
	    {
	      printf("\n\n*** %s does not exist ***\n", reffile);
	      exit(1);
	    }
	}
    }
  
  if ((out_fp = fopen(outfile, filemode))==NULL)
    {
      printf("\n\n*** Error opening %s ***\n", outfile);
      exit(1);
    }

  if (wig && (wig_fp= fopen(wigfile, filemode))==NULL)
    {
      printf("\n\n*** Error opening %s ***\n", wigfile);
      exit(1);
    }
  if (wig && strcmp("",contigfile) && (NT_fp=fopen(contigfile,"r"))==NULL)
    {
      printf("\n\n*** Error opening %s ***\n", contigfile);
      exit(1);
    }

  if(comp && (comp_fp= fopen(compfile, "r"))==NULL)
    {
      printf("\n\n*** Error opening %s ***\n", compfile);
      exit(1);
    }

  if(overlap.print && (overcomp_fp= fopen(compoverlapfile, filemode))==NULL)
    {
      printf("\n\n*** Error opening %s ***\n", wigfile);
      exit(1);
    }

  if(overlap.print && (overreg_fp= fopen(regoverlapfile, filemode))==NULL)
    {
      printf("\n\n*** Error opening %s ***\n", wigfile);
      exit(1);
    }  
}

int platform_strand_endline()
{
  if(!strcmp(filetype,"E"))
    {
      strcpy(pvestring,"F");
      strcpy(negstring,"R");
      return(8);
    }

  else
    {
      
      if(!strcmp(filetype,"M"))
	{return(13);}
      if(!strcmp(filetype,"B"))
	{
	  if(overlap.mode) {return(2);}
	  else             {return(5);}
	}
      if(!strcmp(filetype,"S"))
	{return(9);}
    }
  return(0);
}

struct param get_firstline(int *endfile, char *nextline, struct param par)
{
  if(egets (nextline, 1000, in_fp) == NULL)
    {*endfile = 1;}
  else
    {
      char *nl;
      nl = strchr(nextline, '\n');
      if(nl){*nl = '\0';}
    }

  char c = *nextline;
  char *at = "@";
  char atchar = *at;
  char *hash = "#";
  char hashchar = *hash;
  
  while(c == atchar || c == hashchar)
    {
      if(egets (nextline, 1000, in_fp) == NULL)
	{*endfile = 1;}	  
      else
	{
	  char *nl;
	  nl = strchr(nextline, '\n');
	  if(nl){*nl = '\0';}
	}
      c = *nextline; 
      //printf("%c\n", c);
      //par.totalreads--;
    }

 
  char seqline[1000];
  strcpy(seqline, nextline);
  // ***Calculate read length if unspecified for BED and MAQ*** 
      
  if(par.seqlength == 0 && ((!strcmp(filetype,"B"))||(!strcmp(filetype,"M"))||(!strcmp(filetype,"C"))))
    {
      strcpy(seqline,nextline);
      char *result = NULL;
      result = strtok( seqline, "\t" );
      int i = 0;

      
      while( result != NULL)
	{  
	  //printf(result); printf("\n");
	  if(strcmp(filetype,"M") && i == 1) {par.seqlength = par.seqlength - atoi(result);}
	  if(strcmp(filetype,"M") && i == 2) {par.seqlength = par.seqlength + atoi(result) + 1;}
	  if(!strcmp(filetype,"M") && i == 13) {par.seqlength = atoi(result);}
	  if(!par.pen_inc) {par.pen_inc = par.seqlength;}
	  i++;
	  result = strtok( NULL, "\t" );
	  //printf("%i",par.seqlength);printf("\n");
	}
    }   

  //if(!strcmp(filetype,"S"))  //ignore SAM lines at beginning of file starting with '@'
    //{
      //par.paired = 1;
      //if(par.seqlength == 0)
      //{
      //	  strcpy(seqline,nextline);
      //	  char *result = NULL;
      //  result = strtok( seqline, "\t" );
      //  int i = 0;
      //
      //  while( result != NULL)
      //    { 
      //      if(i == 9) { par.seqlength = strlen(result); }
      //      if(!par.pen_inc) {par.pen_inc = par.seqlength;}
      //      i++;
      //      result = strtok( NULL, "\t" );
      //    }
      //}
  //}

  if(par.fraglength == 0 || par.paired) { par.fraglength = par.seqlength; }
  //if(par.fraglength == 0 && !par.paired) { par.fraglength = par.seqlength; }

  if(!strcmp(filetype,"S"))  //ignore SAM lines at beginning of file starting with '@'
    {
      par.paired = 1;
    }

  if(par.ref)
    {
      if(!par.refseqlength) { par.refseqlength = par.seqlength; }
      if(!par.reffraglength) { par.reffraglength = par.fraglength; }      
    }
  return(par);
}


 void print_output_header(struct param par, char *infile, char *reffile)
{
  //if(par.ref)
  //  {
      fprintf(out_fp, "#Input\t%s\n#Reference\t%s\n#Sequence length\t%d\n#Fragment length\t%d\n#Background\t%f\n#Position Background\t%f\n#Long Background\t%f\n#Threshold\t%f\n#Minimum count above bg\t%d\n#Penalty increase\t%d\n#Quality cutoff\t%f\n#Result cutoff\t%f\n#Penalty factor\t%f\n", infile, reffile, par.seqlength, par.fraglength, par.bg, par.posbg, par.longbg, par.threshold, par.min_above_bg, par.pen_inc, par.qualcutoff, par.resultcutoff, par.refpen);
      fprintf(out_fp, "Region\tStart pos.\tEnd pos.\tCount\tLength\tUnique pos.\tScore\tRef. count\tMax. Coverage\tSummit\n");
      //   }
      // else
      // {
      //    fprintf(out_fp, "#Input\t%s\n#Sequence length\t%d\n#Fragment length\t%d\n#Background\t%f\n#Position Background\t%f\n#Long Background\t%f\n#Threshold\t%d\n#Minimum count above bg\t%d\n#Penalty increase\t%d\n#Quality cutoff\t%f\n#Result cutoff\t%f\n", infile, par.seqlength, par.fraglength, par.bg, par.posbg, par.longbg, par.threshold, par.min_above_bg, par.pen_inc, par.qualcutoff, par.resultcutoff);
      //fprintf(out_fp, "Region\tStart pos.\tEnd pos.\tCount\tLength.\tUnique pos.\tScore\tRef. count\n");
      //}
}



 
void readsampleline(struct readinfo *read,  char *nextline, int *endfile, int fraglength, int seqlength, int qualcutoff, int paired, int overlapmode, int bootstrap) 
{
  int i = 0;
  char *split = NULL;
  char seqline[1000];

  (*read).prevpos = (*read).pos;
  
  //Split first line and assign depending on format
  
  //printf("NEXTLINE:%s\n", nextline);
  split = strtok( nextline, "\t" );

  //printf("PREV:%s\n",(*read).prevchr);
    
  if(!strcmp(filetype,"C"))
    {
      while(split != NULL)
	{
	  switch(i)
	    {       
	    case(0): {strcpy((*read).chr,split); break;}
	    case(1): {(*read).pos = atoi(split); break;}
	    case(2): {(*read).pairlength = atoi(split) - (*read).pos + 1; break;}
	    case(3): {(*read).count = atof(split); break;}
	    case(4): {(*read).pvecount = atof(split); break;}
	    case(5): {(*read).negcount = atof(split); break;}
	    }
	  i++;
	  split = strtok( NULL, "\t" );
	}
      
      if((*read).pvecount)
	{
	  if(!paired) { push((*read).pos,(*read).pos+fraglength-1, (*read).pvecount, &pvepos); }
	  
	  else { push((*read).pos,(*read).pos+(*read).pairlength-1, (*read).pvecount, &pvepos); }
	}
      /*int pc;
	for(pc = 1; pc <= (*read).pvecount; pc++)
	{
	push((*read).pos,(*read).pos+fraglength-1, 1, &pvepos);
	//push(pos + fraglength - 1, &pveend);  //CHANGE
	} */
      
      if(egets (nextline, 1000, in_fp) == NULL)
	{*endfile = 1;}	  
      else
	{
	  char *nl;
	  nl = strchr(nextline, '\n');
	  if(nl){*nl = '\0';}
	}
      
      strcpy(seqline,nextline);
      //printf("SEQ:%s\n", seqline);
      char *nls = NULL;
      nls = strtok( seqline, "\t" );
      if(strcmp((*read).prevchr,(*read).chr)) {(*read).changechr=1;}
    }
  
  else
    {
      (*read).negcount = 0;
  
	  
      for(i=0;i <= endcol; i++)
	{ 
	  if(!strcmp(filetype,"M"))
	    {
	      switch(i)
		{
		case(1): {strcpy((*read).chr,split); break;}
		case(2): {(*read).pos = atoi(split); break;}
		case(3): {strcpy((*read).strand,split); break;}
		case(4): {(*read).pairlength = atoi(split); break;}
		case(5): {(*read).pairflag = atoi(split); break;}
		case(6): {(*read).qual = atoi(split); break;}
		}
	    }
	  
	  else
	    {
	      if(!strcmp(filetype,"B"))
		{
		  switch(i)
		    {
		    case(0): {strcpy((*read).chr,split); break;}
		    case(1): {(*read).pos = atoi(split); break;}
		    case(2): {
		               (*read).pairlength = atoi(split) - (*read).pos + 1; 
			       break;
		             } 
		    case(5): {strcpy((*read).strand,split); break;}
		    }
		}
	      
	      else
		{
		  if(!strcmp(filetype,"E"))
		    {
		      switch(i)
			{
			case(6): {strcpy((*read).chr,split); break;}
			case(7): {(*read).pos = atoi(split); break;}
			case(8): {strcpy((*read).strand,split); break;}
			}
		    }
		  else
		    {
		      if(!strcmp(filetype,"S"))
			{
			  switch(i)
			    {
			    case(1): 
			      {
				if((int)((int)atoi(split) & (int)16))
				  { 
				    strcpy((*read).strand,"-");     
				  }
				else
				  { 
				    strcpy((*read).strand,"+"); 
				  }
				//printf("%s\n", (*read).strand); 780 = 4+8+256+512
				(*read).qual=((int)((int)atoi(split) & (int)780)) ? -1 : 0;
				switch((int)((int)atoi(split) & (int)3))
				  {
				  case(3): { (*read).pairflag = 18; break;}
				  case(1): { (*read).pairflag = 0; break;}
				  default: { (*read).pairflag = -1; break;}
				  }
				break;
			      }   
			    case(2): {strcpy((*read).chr,split); break;}
			    case(3): {(*read).pos = atoi(split); break;}
			    case(4): 
			      {
				if((*read).qual >= 0)
				  {
				    (*read).qual = atoi(split); 
				  }
				if(!strcmp((*read).chr,"*") || (*read).pos == 0)
				  {
				    (*read).qual = -1;
				  }
				break;
			      }
			    case(6): 
			      { 
				if((*read).pairflag == 18 && strcmp(split,"*") && strcmp(split,"="))
				  {
				    (*read).pairflag = 0;
				  }
				break;
			      }
			    case(7):
			      {
				(*read).nextpos = atoi(split);
				if((*read).pairflag == 18 && (*read).nextpos == 0)
				  { 
				    (*read).pairflag = 0; 
				  }
				break;
			      }
			    case(8):
			      {
				(*read).pairlength = atoi(split);
				if((*read).pairflag == 18 && (*read).pairlength <0)
				  {
				    (*read).pairflag = 0;
				  }
				else
				  {
				    if((*read).pairflag == 18 && ((*read).nextpos-(*read).pos) >= (*read).pairlength)
				      {
					(*read).pairflag = 0;
				      }
				  }
				break;
			      }
			    case(9):
			      {
				if((*read).pairflag != 18)
				  {
				    (*read).pairlength = strlen(split);
				  }
				break;
			      }
			    }
			}
		    }
		}
	    }
	  
	  split = strtok( NULL, "\t" );
	}
	
      //printf("ALL1: %s\t%ld\t%s\t%d\t%d\n", (*read).chr, (*read).pos, (*read).strand, (*read).pairlength, (*read).pairflag);
	  
      if(!(*read).pairlength) { (*read).pairlength = seqlength; }

      //if(strcmp(chr,prevchr)) { changechr = 1; }
      struct readinfo nlread;

      strcpy(nlread.chr,(*read).chr);
      nlread.pos = (*read).pos;
      strcpy(nlread.strand,(*read).strand);
      nlread.qual = (*read).qual;
      nlread.pairlength = (*read).pairlength;
      nlread.pairflag = (*read).pairflag;
      int j;
      int firsttime = 1;

	
      //printf("NEXTLINE2:%s\n", nextline);
      if(egets (nextline, 1000, in_fp) == NULL)
	{*endfile = 1;}	  
      else
	{
	  char *nl;
	  nl = strchr(nextline, '\n');
	  if(nl){*nl = '\0';}
	}
      
      // **If quality is high enough(MAQ only), append pvepos if +ve strand
      // **or add 1 to negcount if -ve strand
      
      if(nlread.qual >= qualcutoff)
	{
	  if(!paired) { (*read).negcount += separate(nlread.strand, fraglength, nlread.pos, bootstrap); }
	  //else{ if(nlread.pairflag  == 18) { (*read).negcount += separate(nlread.strand, nlread.pairlength, nlread.pos); }}
	  else{ 
	    if(nlread.pairflag == -1) {
	      (*read).negcount += separate(nlread.strand, (fraglength > nlread.pairlength) ? fraglength : nlread.pairlength, nlread.pos, bootstrap); }
	    if(nlread.pairflag == 18) { 
	      separate(nlread.strand, nlread.pairlength, nlread.pos, bootstrap); } }
	      //if(fraglength){ separate(nlread.strand, fraglength, nlread.pos, bootstrap); }
	      //else { separate(nlread.strand, nlread.pairlength, nlread.pos, bootstrap); } } }
	  // ***If paired ignore negative strand reads (except for count file format)***
	}
      


      // **Continue reading in file while position and chromosome are the same
      while(nlread.pos == (*read).pos && !strcmp(nlread.chr,(*read).chr) && !(*endfile))
	{ 	 
	  if(!firsttime)
	    {
	      if(!(*endfile) && egets (nextline, 1000, in_fp) == NULL)
		{*endfile = 1; break;}
	      else
		{
		  char *nl;
		  nl = strchr(nextline, '\n');
		  if(nl){*nl = '\0';}
		}
	      if(nlread.qual >= qualcutoff)
		{
		  if(!paired) { (*read).negcount += separate(nlread.strand, fraglength, nlread.pos, bootstrap); }
		  else{
		    if(nlread.pairflag == -1) {
		      (*read).negcount += separate(nlread.strand, (fraglength > nlread.pairlength) ? fraglength : nlread.pairlength, nlread.pos, bootstrap); }
		    if(nlread.pairflag  == 18) { separate(nlread.strand, nlread.pairlength, nlread.pos, bootstrap); } }
		  // ***If paired ignore negative strand reads (except for count file format)***
		}
	    }
	  firsttime = 0;
	  
	  
	  //printf("NEXTLINE3:%s\n", nextline);
	  strcpy(seqline,nextline);
	  char *nls = NULL;
	  nls = strtok( seqline, "\t" );
	  
	  //  if(nlsqual >= qualcutoff)
	  //{ separate(nlsstrand);}
	  
	  for(j = 0; j<=endcol; j++)
	    {
	      if(!strcmp(filetype,"M"))
		{
		  switch(j)
		    {
		    case(1): {strcpy(nlread.chr,nls); break;}
		    case(2): {nlread.pos = atoi(nls); break;}
		    case(3): {strcpy(nlread.strand,nls); break;}
		    case(4): {nlread.pairlength = atoi(nls); break;}
		    case(5): {nlread.pairflag = atoi(nls); break;}
		    case(6): {nlread.qual = atoi(nls); break;}
		    }
		}
		
	      else
		{
		  if(!strcmp(filetype,"B"))
		    {
		      switch(j)
			{
			case(0): {strcpy(nlread.chr,nls); break;}
			case(1): {nlread.pos = atoi(nls); break;}
			case(2): {
			          nlread.pairlength = atoi(nls) - nlread.pos + 1;
				  if(overlapmode) { i=5;}
				  break;
			         }  
			case(5): {strcpy(nlread.strand,nls); break;}
			}
		    }
		  
		  else
		    {
		      if(!strcmp(filetype,"E"))
			{
			  switch(j)
			    {
			    case(6): {strcpy(nlread.chr,nls); break;}
			    case(7): {nlread.pos = atoi(nls); break;}
			    case(8): {strcpy(nlread.strand,nls); break;}
			    }
			}
		      else
			{
			  if(!strcmp(filetype,"S"))
			    {
			      switch(j)
				{
				
				case(1): 
				  {
				    if((int)((int)atoi(nls) & (int)16))
				      { 
					strcpy(nlread.strand,"-");     
				      }
				    else
				      { 
					strcpy(nlread.strand,"+"); 
				      }
				    //printf("%s\n", nlread.strand); 780 = 4+8+256+512
				    nlread.qual=((int)((int)atoi(nls) & (int)780)) ? -1 : 0;
				    switch((int)((int)atoi(nls) & (int)3))
				      {
				      case(3): { nlread.pairflag = 18; break;}
				      case(1): { nlread.pairflag = 0; break;}
				      default: { nlread.pairflag = -1; break;}
				      }
				    break;
				  }   
				case(2): {strcpy(nlread.chr,nls); break;}
				case(3): {nlread.pos = atoi(nls); break;}
				case(4): 
				  {
				    if(nlread.qual >= 0)
				      {
					nlread.qual = atoi(nls);
				      }
				    if(!strcmp(nlread.chr,"*") || nlread.pos == 0)
				      {
					nlread.qual = -1;
				      }
				    break;
				  }
				case(6): 
				  { 
				    if(nlread.pairflag == 18 && strcmp(nls,"*") && strcmp(nls,"="))
				      {
					nlread.pairflag = 0;
				      }
				    break;
				  }
				case(7):
				  {
				    nlread.nextpos = atoi(nls);
				    if(nlread.pairflag == 18 && nlread.nextpos == 0)
				      { 
					nlread.pairflag = 0; 
				      }
				    break;
				  }
				case(8):
				  {
				    nlread.pairlength = atoi(nls);
				    if(nlread.pairflag == 18 && nlread.pairlength <0)
				      {
					nlread.pairflag = 0;
				      }
				    else
				      {
					if(nlread.pairflag == 18 && (nlread.nextpos-nlread.pos) >= nlread.pairlength)
					  {
					    nlread.pairflag = 0;
					  }
				      }
				    break;
				  }
				case(9):
				  {
				    if(nlread.pairflag != 18)
				      {
					nlread.pairlength = strlen(nls);
				      }
				    break;
				  }
				}
			    }
			}
		    }
		}
		
	      nls = strtok( NULL, "\t" );
	    }
	  //printf("NEXTLINE4:%s\n", nextline);	    
	  //printf("ALL: %s\t%ld\t%s\t%d\t%d\n", nlread.chr, nlread.pos, nlread.strand, nlread.pairlength, nlread.pairflag);
	  
	  //if(!endfile && egets (nextline, 1000, in_fp) == NULL)
	  //{endfile = 1;}
	  
	}
	  
      //printf("B1%s\t%s\t%s\n",(*read).prevchr,(*read).chr, nlread.chr);
      if(strcmp(nlread.chr,(*read).chr)) {(*read).changechr=1;}    
    }
  //return(read);
}


void print_out_lastregion(struct region *lastregion, double binom_p, int peakversion)
{
  
  //printf("%s\t%ld\t%ld\t%ld\t%d\t%f\t%f\t%ld\n", (*lastregion).chr, (*lastregion).start, (*lastregion).end, (*lastregion).end-(*lastregion).start+1, (*lastregion).pos, (*lastregion).print_total, (*lastregion).peakheightmax, (long int)floor((*lastregion).summit));
  char halftext[3] = "";
  if((double)floor((*lastregion).summit) != (*lastregion).summit)
    { strcpy(halftext, ".5"); }

  if(peakversion)
    {
      (*lastregion).start = (*lastregion).seen_start;
      (*lastregion).end = (*lastregion).seen_end;
    }

  // ***Unless using countfile (double) converted counts to integers***
  if(strcmp(filetype,"C"))  
    {
      int intcount = floor((*lastregion).count + 0.5);
      int intrefcount = floor((*lastregion).refcount + 0.5);
      fprintf(out_fp,"%s\t%ld\t%ld\t%d\t%ld\t%d\t%f\t%d\t%f\t%ld%s\t%f\n", (*lastregion).chr, (*lastregion).start, (*lastregion).end, intcount, (*lastregion).end-(*lastregion).start+1, (*lastregion).pos, (*lastregion).print_total, intrefcount, (*lastregion).peakheightmax, (long int)floor((*lastregion).summit), halftext, 1-pbinom(intcount,intcount+intrefcount,binom_p));   
      //printf("%s\t%ld\t%ld\t%d\t%ld\t%d\t%f\t%d\t%f\t%ld%s\n", (*lastregion).chr, (*lastregion).start, (*lastregion).end, intcount, (*lastregion).end-(*lastregion).start+1, (*lastregion).pos, (*lastregion).print_total, intrefcount, (*lastregion).peakheightmax, (long int)floor((*lastregion).summit), halftext);    
    }

  else
    {      
      fprintf(out_fp,"%s\t%ld\t%ld\t%f\t%ld\t%d\t%f\t%f\t%f\t%ld%s\n", (*lastregion).chr, (*lastregion).start, (*lastregion).end, (*lastregion).count, (*lastregion).end-(*lastregion).start+1, (*lastregion).pos, (*lastregion).print_total,(*lastregion).refcount, (*lastregion).peakheightmax, (long int)floor((*lastregion).summit), halftext);
    }
  //fprintf(out_fp,"%s\t%ld\t%ld\t%f\t%ld\t%d\t%f\t%f\n", (*lastregion).chr, (*lastregion).start, (*lastregion).end, (*lastregion).count, (*lastregion).end-(*lastregion).start+1, (*lastregion).pos, (*lastregion).total,(*lastregion).refcount);
  //printf("%s\t%ld\t%ld\t%f\t%ld\t%d\t%f\t%f\n", (*lastregion).chr, (*lastregion).start, (*lastregion).end, (*lastregion).count, (*lastregion).end-(*lastregion).start+1, (*lastregion).pos, (*lastregion).total, (*lastregion).refcount);
}



//***If zipped use zlib open and read functions***
FILE *eopen(const char *path, const char *mode)
{
  if(zip)
    {
      return gzopen(path, mode);
    }
  else
    {
      return fopen(path, mode);
    }
}

char *egets(char *line, int n, FILE *fp)
{
  if(zip)
    {
      return gzgets(fp, line, n);
    }
  else
    {
      return fgets(line, n, fp);
    }
}

long int count_sample_lines(struct param *par, char infile[1000], int samplecount)
{
  FILE *count_fp;
  char countline[1000];     
  long int nqreads = 0;
  char *at = "@";
  char atchar = *at;
  
  char *hash = "#";
  char hashchar = *hash;
	  

  //if(!strcmp(filetype,"F"))
  //  {
  //    char countcommand[1100] = "samtools flagstat ";
  //    strcat(countcommand, infile);
  //    strcat(countcommand, " | awk '{if(NR==1){print $1}}'");
      
  //    count_fp = popen(countcommand, "r");
      
  //    char countline[1000];
  //    fgets(countline, 1000, count_fp);
  //    char *nl;
  //    nl = strchr(countline, '\n');
  //    if(nl){*nl = '\0';}
  //    
  //    nqreads = atoi(countline);
  //  }


  if(zip)
    { count_fp = gzopen(infile, "r"); }
  else
    { 
      if(strcmp(filetype,"F"))
	{
	  count_fp = fopen(infile,"r"); 
	}
      else
	{
	  char samtoolscommand[1000] = "samtools view ";
	  strcat(samtoolscommand, infile);
	  count_fp = popen(samtoolscommand, "r");
	}
    }
      

     
  if(!strcmp(filetype, "F") || !strcmp(filetype,"S"))
    {
      while(egets (countline, 1000, count_fp) != NULL)
	{ 
	  struct readinfo countread;
	  char countchar = *countline;
	  if(countchar == atchar){ continue; }
	  
	  char *split = NULL;
	  split = strtok( countline, "\t" );
	    
	  split = strtok( NULL, "\t");

	  if((int)((int)atoi(split) & (int)16))
	    { 
	      strcpy(countread.strand,"-");     
	    }
	  else
	    { 
	      strcpy(countread.strand,"+"); 
	    }
	  //printf("%s\n", countread.strand); 780 = 4+8+256+512
	  countread.qual=((int)((int)atoi(split) & (int)780)) ? -1 : 0;

	  if((int)((int)atoi(split) & (int)3) == 1)
	    {
	      continue;
	    }
	  countread.pairflag = ((int)((int)atoi(split) & (int)3) == 3) ? 18 : -1;
	      	  

	  split = strtok( NULL, "\t");
	  strcpy(countread.chr,split); 
	  
	  split = strtok( NULL, "\t");
	  countread.pos = atoi(split);

	  split = strtok( NULL, "\t");
	  if(countread.qual >= 0)
	    {countread.qual = atoi(split);}

	  if(countread.qual == 0 || !strcmp(countread.chr,"*") || countread.pos == 0)
	    {
	      continue;
	    }

	  split = strtok( NULL, "\t");
	  split = strtok( NULL, "\t");
	  
	  if(countread.pairflag == 18 && strcmp(split,"*") && strcmp(split,"="))
	    {
	      continue;
	    }

	  split = strtok( NULL, "\t");
	  countread.nextpos = atoi(split);
	  if(countread.pairflag == 18 && countread.nextpos == 0)
	    { 
	      continue; 
	    }

	  split = strtok( NULL, "\t");
	  countread.pairlength = atoi(split);
	  if(countread.pairflag == 18 && (countread.pairlength <0 || (countread.nextpos-countread.pos) >= countread.pairlength))
	    {
	      continue;
	    }
	  
	  if(countread.pairflag != 0)
	    {
	      if(countread.qual >= (*par).qualcutoff)
		{ nqreads++;}
	      if(samplecount)
		{ (*par).nlines++; }
	    }
	}
    }
  
  else
    {
      if(strcmp(filetype,"M") || ((*par).qualcutoff == 0 && (*par).paired == 0))
	{
	  
	  while(egets (countline, 1000, count_fp) != NULL)
	    { 
	      char countchar = *countline;
	      //printf("%c %c\n", countchar, atchar);
	      if(countchar != atchar && countchar != hashchar)
		{ nqreads++; }
	    }
	  if(samplecount)
	    { (*par).nlines = nqreads; }
	}
    
      
      else
	{  
	  char *c = NULL;
	  int i;
	  
	  if((*par).paired)
	    {
	      while(egets (countline, 1000, count_fp) != NULL)
		{
		  i = 1;
		  c = strtok(countline, "\t");
		  while( i++ <= 3)
		    { c = strtok(NULL, "\t"); }
		  //printf("s %s\n",c);
		  if(strcmp(c,"+"))
		    { 
		      //printf( "Neg\n"); 
		      continue; 
		    }
		  while( i++ <= 6)
		    { c = strtok(NULL, "\t"); }
		  if(atoi(c) != 18)  //Pair flag
		    { 
		      //printf ("Flag\n"); 
		      continue; 
		    }
		  c = strtok(NULL, "\t");
		  if(atoi(c) >= (*par).qualcutoff)
		    { 
		      //printf ("Qual\n"); 
		      nqreads++; 
		    }
		  if(samplecount)
		    { (*par).nlines++; }
		}
	    }
	  else
	    {	  
	      while(egets (countline, 1000, count_fp) != NULL)
		{	      
		  i = 1;
		  c = strtok(countline, "\t");
		  while( i++ <= 7)
		    { c = strtok(NULL, "\t"); }
		  //printf("s %s\n",c);
		  if(atoi(c) >= (*par).qualcutoff)
		    { nqreads++;}
		  if(samplecount)
		    { (*par).nlines++; }		
		}
	    }
	}	  
    }
  if(zip)
    { gzclose(count_fp); }
  else
    { fclose(count_fp); }
  //printf("%ld total reads.\n", nqreads);
  return nqreads;
}





void readrefline(struct readinfo *refread, int *refendfile, int reffraglength, int paired, int bootstrap)
{
  strcpy((*refread).prevchr, (*refread).chr);
  char refline[1000];

  if(egets (refline, 1000, ref_fp) == NULL)
    { *refendfile = 1; }
  else
    {
      char *refnl;
      refnl = strchr(refline, '\n');
      if(refnl){*refnl = '\0';}
    }

  //if(!strcmp(filetype,"S"))
  // {
      char c = *refline;
      
      char *at = "@";
      char atchar = *at;

      char *hash = "#";
      char hashchar = *hash;

      while(c == atchar || c == hashchar)
	{
	  if(egets (refline, 1000, in_fp) == NULL)
	    {*refendfile = 1;}	  
	  else
	    {
	      char *nl;
	      nl = strchr(refline, '\n');
	      if(nl){*nl = '\0';}
	    }
	  c = *refline; 
	  //par.totalrefreads--;
	}

      //}
  

  //printf("REF: %s\n",refline);

  if(!(*refendfile))
    {
      int ri = 0;
      char *refsplit = NULL;

      refsplit = strtok( refline, "\t" );

      for(ri=0;ri <= endcol; ri++)
	{ 
	  
	  if(!strcmp(filetype,"C"))
	    {
	      switch(ri)
		{       
		case(0): {strcpy((*refread).chr,refsplit); break;}
		case(1): {(*refread).pos = atoi(refsplit); break;}
		case(2): {(*refread).pairlength = atoi(refsplit) - (*refread).pos + 1; break;}
		case(3): {(*refread).count = atoi(refsplit); break;}
		case(4): {(*refread).pvecount = atoi(refsplit); break;}
		case(5): {(*refread).negcount = atoi(refsplit); break;}
		}
	    }
	  
	  
	  else
	    {
	      if(!strcmp(filetype,"M"))
		{
		  switch(ri)
		    {
		    case(1): {strcpy((*refread).chr,refsplit); break;}
		    case(2): {(*refread).pos = atoi(refsplit); break;}
		    case(3): {strcpy((*refread).strand,refsplit); break;}
		    case(4): {(*refread).pairlength = atoi(refsplit); break;}
		    case(5): {(*refread).pairflag = atoi(refsplit); break;}
		    case(6): {(*refread).qual = atoi(refsplit); break;}
		    }
		}
		  
	      else
		{
		  if(!strcmp(filetype,"B"))
		    {
		      switch(ri)
			{
			case(0): {strcpy((*refread).chr,refsplit); break;}
			case(1): {(*refread).pos = atoi(refsplit); break;}
			case(2): {(*refread).pairlength = atoi(refsplit) - (*refread).pos + 1; break;}
			case(5): {strcpy((*refread).strand,refsplit); break;}
			}
		    }

		  else
		    {
		      if(!strcmp(filetype,"E"))
			{
			  switch(ri)
			    {
			    case(6): {strcpy((*refread).chr,refsplit); break;}
			    case(7): {(*refread).pos = atoi(refsplit); break;}
			    case(8): {strcpy((*refread).strand,refsplit); break;}
			    }
			}
		      else
			{
			  if(!strcmp(filetype,"S"))
			    {
			      switch(ri)
				{
				case(1): 
				  {
				    if((int)((int)atoi(refsplit) & (int)16))
				      { 
					strcpy((*refread).strand,"-");     
				      }
				    else
				      { 
					strcpy((*refread).strand,"+"); 
				      }
				    //printf("%s\n", (*refread).strand); 780 = 4+8+256+512
				    (*refread).qual=((int)((int)atoi(refsplit) & (int)780)) ? -1 : 0;
				    switch((int)((int)atoi(refsplit) & (int)3))
				      {
				      case(3): { (*refread).pairflag = 18; break;}
				      case(1): { (*refread).pairflag = 0; break;}
				      default: { (*refread).pairflag = -1; break;}
				      }
				    break;
				  }   
				case(2): {strcpy((*refread).chr,refsplit); break;}
				case(3): {(*refread).pos = atoi(refsplit); break;}
				case(4): 
				  {
				    if((*refread).qual >= 0)
				      {
					(*refread).qual = atoi(refsplit);  
				      }
				    if(!strcmp((*refread).chr,"*") || (*refread).pos == 0)
				      {
					(*refread).qual = -1;
				      }
				    break;
				  }
				case(6): 
				  { 
				    if((*refread).pairflag == 18 && strcmp(refsplit,"*") && strcmp(refsplit,"="))
				      {
					(*refread).pairflag = 0;
				      }
				    break;
				  }
				case(7):
				  {
				    (*refread).nextpos = atoi(refsplit);
				    if((*refread).pairflag == 18 && (*refread).nextpos == 0)
				      { 
					(*refread).pairflag = 0; 
				      }
				    break;
				  }
				case(8):
				  {
				    (*refread).pairlength = atoi(refsplit);
				    if((*refread).pairflag == 18 && (*refread).pairlength <0)
				      {
					(*refread).pairflag = 0;
				      }
				    else
				      {
					if((*refread).pairflag == 18 && ((*refread).nextpos-(*refread).pos) >= (*refread).pairlength)
					  {
					    (*refread).pairflag = 0;
					  }
				      }
				    break;
				  }
				case(9):
				  {
				    if((*refread).pairflag != 18)
				      {
					(*refread).pairlength = strlen(refsplit);
				      }
				    break;
				  }
				}
			    }
			}
		    }
		}		      
	    }
	  refsplit = strtok( NULL , "\t");
	}
    }
  //printf("RALL1: %s\t%ld\t%s\t%d\t%d\n", (*refread).chr, (*refread).pos, (*refread).strand, (*refread).pairlength, (*refread).pairflag);
  if(!paired){ (*refread).pairlength = reffraglength; }
  if(strcmp(filetype,"C"))
    {
      int poi = rpois(bootstrap,1);
      //(*refread).count = 1;
      (*refread).count = poi;
      if(!strcmp((*refread).strand,negstring))
	{ (*refread).pvecount = 0; (*refread).negcount = poi; }
	//{ (*refread).pvecount = 0; (*refread).negcount = 1; }
      else
	{
	  if(!strcmp((*refread).strand,pvestring))
	    { (*refread).pvecount = poi; (*refread).negcount = 0; }
	    //{ (*refread).pvecount = 1; (*refread).negcount = 0; }
	  else{printf("Reference strand not recognised.\n");exit(1);}
	}

      if(paired && ((!strcmp((*refread).strand,negstring)) || (*refread).pairflag != 18))
	{ (*refread).count = 0; (*refread).pvecount = 0; (*refread).negcount = 0; }
    }
  //return(ref);
}
