/*  Copyright 2013 Daniel Wilson and Xavier Didelot.
 *
 *  xmfa.h
 *  Part of ClonalFrameML
 *
 *  The myutils library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  The myutils library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with the myutils library. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#define UNLINKED 9
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

char
convert (char in)
{
  switch (in)
    {
    case 'a':
      return '0';
    case 'A':
      return '0';
    case 't':
      return '1';
    case 'T':
      return '1';
    case 'c':
      return '2';
    case 'C':
      return '2';
    case 'g':
      return '3';
    case 'G':
      return '3';
    default:
      return 'N';
    }
}

void
Fgets (char *buf, int l, FILE * f)
{
  fgets (buf, l, f);
  while (strncmp (buf, "#", 1) == 0)
    fgets (buf, l, f);
  if (buf[strlen (buf) - 1] == '\n')
    buf[strlen (buf) - 1] = '\0';
  if (buf[strlen (buf) - 1] == '\r')
    buf[strlen (buf) - 1] = '\0';
}

void
extractName (char *buf, char *buf2)
{
  if (buf[0] != '>')
    {
      buf2[0] = '\0';
      return;
    }
  int i = 1;
  while (buf[i] == ' ')
    i++;
  strcpy (buf2, buf + i);
  if (strstr (buf2, ":") != NULL)
    (strstr (buf2, ":"))[0] = '\0';
}


DNA readXMFA(const char *filename) {
  printf ("Reading xmfa file...\n");
  char *names[10000];
  char buf[10000];
  char buf2[1000];
  //First count the number of genomes
  FILE *f = fopen (filename, "r");
  if (f == NULL)
    {
      printf ("Unable to open input file\n");
      abort ();
    }
  int nbgenomes = 0;
  int i, l, j;
  int deb = 0;
  while (!feof (f))
    {
      Fgets (buf, 10000, f);
      if (buf[0] == '>')
	{
	  names[nbgenomes++] = (char *) calloc (1000, sizeof (char));
	  extractName (buf, names[nbgenomes - 1]);
	}
      if (buf[0] == '=' || feof (f))
	break;
    }

  //Then count the number of LCBs where all sequences are present
  rewind (f);
  int *check = (int *) calloc (nbgenomes, sizeof (int));
  int lcb = 0;
  while (!feof (f))
    {
      Fgets (buf, 10000, f);
      if (buf[0] == '>')
	{
	  extractName (buf, buf2);
	  for (int i = 0; i < nbgenomes; i++)
	    if (strcmp (buf2, names[i]) == 0)
	      {
		check[i] = 1;
		break;
	      }
	}
      if (strncmp (buf, "=", 1) == 0 || feof (f))
	{
	  int sum = 0;
	  for (int i = 0; i < nbgenomes; i++)
	    sum += check[i];
	  if (sum == nbgenomes)
	    {
	      lcb++;
	      for (int i = 0; i < nbgenomes; i++)
		check[i] = 0;
	    }
	  else
	    {
	      int sum = 0;
	      for (int i = 0; i < nbgenomes; i++)
		sum += check[i];
	      if (sum > 0)
		{
		  int missing = 0;
		  for (int i = 0; i < nbgenomes; i++)
		    if (check[i] == 0)
		      {
			missing = i;
			break;
		      }
		  printf
		    ("Warning: fragment %d is incomplete (does not contain strain %s)\nWarning: ClonalFrame is running only on the first %d fragments\n",
		     lcb + 1, names[missing], lcb);
		}
	      break;
	    }
	}
    }
  free (check);
  //Then count the length of the alignment
  rewind (f);
  int length = 0;
  for (l = 0; l < lcb; l++)
    {
      while (1)
	{
	  Fgets (buf, 10000, f);
	  if (buf[0] == '>')
	    {
	      extractName (buf, buf2);
	      if (strcmp (buf2, names[0]) == 0)
		break;
	    }
	}
      Fgets (buf, 10000, f);
      while (buf[0] != '>' && buf[0] != '=' && !feof (f))
	{
	  length += strlen (buf);
	  Fgets (buf, 10000, f);
	}
    }
  printf ("Found %d strains, %d blocks and %d sites.\n", nbgenomes, lcb,
	  length);
  fclose (f);

  //Then read the data
  char **data = (char **) calloc (nbgenomes, sizeof (char *));
  for (int i = 0; i < nbgenomes; i++)
    data[i] = (char *) calloc (length + lcb - 1, sizeof (char));
  FILE **in = (FILE **) calloc (nbgenomes, sizeof (FILE *));
  char **bufin = (char **) calloc (nbgenomes, sizeof (char *));
  int *pos = (int *) calloc (nbgenomes, sizeof (int));
  for (i = 0; i < nbgenomes; i++)
    {
      in[i] = fopen (filename, "r");
      bufin[i] = (char *) calloc (10000, sizeof (char));
    }
  for (l = 0; l < lcb; l++)
    {				//For each LCB
      for (i = 0; i < nbgenomes; i++)
	{			//For each genome
	  //Find next bit of alignment
	  extractName (bufin[i], buf2);
	  while (strcmp (names[i], buf2) != 0)
	    {
	      Fgets (bufin[i], 10000, in[i]);
	      extractName (bufin[i], buf2);
	    }
	  if (i == 0)
	    {
	      if (strstr (bufin[i], ":") != NULL)
		sscanf (strstr (bufin[i], ":") + 1, "%d", &deb);
	      else
		deb++;
	    }
	  Fgets (bufin[i], 10000, in[i]);
	  //And copy it
	  while (bufin[i][0] != '>' && bufin[i][0] != '=' && !feof (in[i]))
	    {
	      for (j = 0; j < strlen (bufin[i]); j++)
		data[i][pos[i]++] = convert (bufin[i][j]);
	      Fgets (bufin[i], 10000, in[i]);
	    }
	  //Add the unlinked symbol if there are more LCB coming
	  if (l < lcb - 1)
	    data[i][pos[i]++] = UNLINKED;
	}
    }
  for (i = 0; i < nbgenomes; i++)
    if (pos[i] != length + lcb - 1)
      {
	printf ("Invalid input file: the sequence for %s is incomplete.\n",
		names[i]);
	abort ();
      }

  //Clear the memory
  for (i = 0; i < nbgenomes; i++)
    {
      fclose (in[i]);
      free (bufin[i]);
    }
  free (in);
  free (pos);
  free (bufin);

DNA fa;

//Clear memory
  for (i = 0; i < nbgenomes; i++)
    free (names[i]);
  for (i = 0; i < nbgenomes; i++)
    free (data[i]);
  free (data);
printf("All done.\n");
return fa;
}

