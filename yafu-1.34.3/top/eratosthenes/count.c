/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

       				   --bbuhrow@gmail.com 7/1/10
----------------------------------------------------------------------*/

#include "soe.h"

uint64 count_line(soe_staticdata_t *sdata, uint32 current_line)
{
	//extract stuff from the thread data structure
	uint8 *line = sdata->lines[current_line];
	uint64 numlinebytes = sdata->numlinebytes;
	uint64 lowlimit = sdata->lowlimit;
	uint64 prodN = sdata->prodN;
	uint64 *flagblock64 = (uint64 *)line;
	uint8 *flagblock;
	uint64 i, k, it = 0;
	int ix;
	int done, kx;
	uint64 prime;

	for (i=0;i<(numlinebytes >> 3);i++)
	{
		/* Convert to 64-bit unsigned integer */    
		uint64 x = flagblock64[i];
		    
		/*  Employ bit population counter algorithm from Henry S. Warren's
			*  "Hacker's Delight" book, chapter 5.   Added one more shift-n-add
			*  to accomdate 64 bit values.
			*/
		x = x - ((x >> 1) & 0x5555555555555555ULL);
		x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
		x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
		x = x + (x >> 8);
		x = x + (x >> 16);
		x = x + (x >> 32);

		it += (x & 0x000000000000003FULL);
	}

	//potentially misses the last few bytes
	//use the simpler baseline method to get these few
	flagblock = line;
	for (k=0; k<(numlinebytes & 7);k++)
	{
		it += (flagblock[(i<<3)+k] & nmasks[0]) >> 7;
		it += (flagblock[(i<<3)+k] & nmasks[1]) >> 6;
		it += (flagblock[(i<<3)+k] & nmasks[2]) >> 5;
		it += (flagblock[(i<<3)+k] & nmasks[3]) >> 4;
		it += (flagblock[(i<<3)+k] & nmasks[4]) >> 3;
		it += (flagblock[(i<<3)+k] & nmasks[5]) >> 2;
		it += (flagblock[(i<<3)+k] & nmasks[6]) >> 1;
		it += (flagblock[(i<<3)+k] & nmasks[7]);
	}

	//eliminate the primes flaged that are above or below the
	//actual requested limits. as both of these can change to 
	//facilitate sieving we'll need to compute them, and
	//decrement the counter if so.
	//this is a scarily nested loop, but it only should iterate
	//a few times.
	done = 0;
	for (ix=(int)numlinebytes-1;ix>=0 && !done;ix--)
	{
		for (kx=BITSINBYTE-1;kx>=0;kx--)
		{
			if (line[ix] & nmasks[kx])
			{
				prime = prodN * ((uint64)ix * (uint64)BITSINBYTE + (uint64)kx) + 
					(uint64)sdata->rclass[current_line] + lowlimit;
				if (prime > sdata->orig_hlimit)
					it--;
				else
				{
					done = 1;
					break;
				}
			}
		}
	}
	done = 0;
	for (ix=0;ix<numlinebytes && !done;ix++)
	{
		for (kx=0;kx<8;kx++)
		{
			if (line[ix] & nmasks[kx])
			{
				prime = prodN * ((uint64)ix * (uint64)BITSINBYTE + (uint64)kx) + 
					(uint64)sdata->rclass[current_line] + lowlimit;
				if (prime < sdata->orig_llimit)
					it--;
				else
				{
					done = 1;
					break;
				}
			}
		}
	}

	return it;
}

void count_line_special(thread_soedata_t *thread_data)
{
	//extract stuff from the thread data structure
	soe_staticdata_t *sdata = &thread_data->sdata;
	uint32 current_line = thread_data->current_line;
	uint8 *line = sdata->lines[current_line];	
	uint64 numlinebytes = sdata->numlinebytes;
	uint64 lowlimit = sdata->lowlimit;
	uint64 prodN = sdata->prodN;
	uint64 *flagblock64 = (uint64 *)line;
	uint64 i, k, it, lower, upper;
	int ix;
	int64 start, stop;

	//zero out any bits below the requested range
	for (i=lowlimit + sdata->rclass[current_line], ix=0; i < sdata->orig_llimit; i += prodN, ix++)
		line[ix >> 3] &= masks[ix & 7];
	
	//and any high bits above the requested range
	for (i=sdata->highlimit + sdata->rclass[current_line] - prodN, ix=0; i > sdata->orig_hlimit; i -= prodN, ix++)
		line[numlinebytes - 1 - (ix >> 3)] &= masks[7 - (ix & 7)];

	//count each block of 1e9
	lower = sdata->orig_llimit;
	upper = lower;
	k = 0;
	thread_data->linecount = 0;
	while (upper != sdata->orig_hlimit)
	{
		//set the bounds for the next batch
		upper = upper + 1000000000; 
		if (upper > sdata->orig_hlimit)
			upper = sdata->orig_hlimit;

		//find the starting byte number.  first find the number of bits between the current lower
		//limit and the start of the line.
		start = (int64)((lower - lowlimit) / prodN);

		//we'll be counting in 64 bit chunks, so compute how many 64 bit chunks this is
		start /= 64;
		
		//start a little before the range, to account for rounding errors
		start -= 2;

		if (start < 0) start = 0;

		//find the stopping byte number: first find the number of bits between the current upper
		//limit and the start of the line.
		stop = (int64)((upper - lowlimit) / prodN);

		//we'll be counting in 64 bit chunks, so compute how many 64 bit chunks this is
		stop /= 64;
		
		//stop a little after the range, to account for rounding errors
		stop += 2;

		if (stop > (int64)(numlinebytes >> 3)) 
			stop = (numlinebytes >> 3);

		//count these bytes
		it = 0;
		for (int64 ix = start; ix < stop; ix++)
		{
			/* Convert to 64-bit unsigned integer */    
			uint64 x = flagblock64[ix];
		    
			/*  Employ bit population counter algorithm from Henry S. Warren's
			 *  "Hacker's Delight" book, chapter 5.   Added one more shift-n-add
			 *  to accomdate 64 bit values.
			 */
			x = x - ((x >> 1) & 0x5555555555555555ULL);
			x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
			x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
			x = x + (x >> 8);
			x = x + (x >> 16);
			x = x + (x >> 32);

			it += (x & 0x000000000000003FULL);
		}

		//then correct the counts
		//zero out any bits below the requested range
		for (i= (start * 64) * prodN + sdata->rclass[current_line] + lowlimit, ix=0; i < lower; i += prodN, ix++)
		{
			if (line[(ix >> 3) + (start << 3)] & nmasks[ix & 7])
				it--;
		}
		
		//and any high bits above the requested range
		for (i=(stop * 64) * prodN + sdata->rclass[current_line] - prodN + lowlimit, ix=0; i > upper; i -= prodN, ix++)
		{
			if (line[(stop << 3) - 1 - (ix >> 3)] & nmasks[7 - (ix & 7)])
				it--;

		}
		
		//add the count to the special array
		thread_data->ddata.special_count[k] = (uint32)it;
		thread_data->linecount += it;
		k++;
		lower = upper;
	}

}
