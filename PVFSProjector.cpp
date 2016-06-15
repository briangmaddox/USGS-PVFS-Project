//Implementation file for the PVFSProjector
//By Chris Bilderback

#include "PVFSProjector.h"

//********************************************************************
PVFSProjector::PVFSProjector() : MpiProjector()
{}

//********************************************************************
PVFSProjector::~PVFSProjector()
{}

//********************************************************************
void PVFSProjector::setPartitionNumber(const unsigned int & inPart) throw()
{
  //resize the vectors
  mcounters.resize(inPart);
  mstop.resize(inPart);
}

//********************************************************************
unsigned int PVFSProjector::getPartitionNumber() const throw()
{
  return mcounters.size();
}


//********************************************************************
void PVFSProjector::writeToPVFS(bool inWriteRaw) throw()
{
  //use slave local for this and assume that they are useing the 
  //correct slave that will actually  write to pvfs.
  slavelocal = inWriteRaw;
}


//********************************************************************
void PVFSProjector::project(BaseProgress * progress = NULL) 
  throw(ProjectorException)
{
  int mytid(0);                                  //my pvm id
  PmeshLib::ProjectionMesh * pmesh(0);           //projection mesh
  int counter(0);
  int scounter(0);

  try
  {
    if (!mcounters.size())           //must have a partition num  
    {
      Projector::project(progress);
      return;
    }

    if (static_cast<unsigned int>(numofslaves) < mcounters.size()) 
    {
      Projector::project(progress);
      return;
    }

    
    if (!fromprojection || !toprojection)        //check for projection
    {
      throw ProjectorException(PROJECTOR_ERROR_UNKOWN);
    }
    
    pmesh = setupForwardPmesh();                 //try setup the forward
                                                 //pmesh
    
    getExtents(pmesh);                           //get the extents
    
    
    if(cache)                                    //delete the cache
    {
      delete cache;
      cache = NULL;
    }
   
    if(!slavelocal)
    {
      setupOutput(outfile);                        //create the output file
    }
    else
    {
      setUpRawOutput();
    }
    
    if (pmesh)                                     //delete uneeded mesh
    {
      delete pmesh;
      pmesh = NULL;
    }

    //check the rank in MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mytid);

    if (mytid != 0)
    {
      //master must be rank zero
      throw ProjectorException(PROJECTOR_ERROR_UNKOWN);
    }
    
    
    //figure out the mapping
    scounter = -1;
    for(counter = 1; counter <= numofslaves; ++counter)
    {
      if (counter % numofslaves/mcounters.size() == 0)
        ++scounter;
      
      membership[counter] = scounter;
    }

    //figure out the start stop stuff
    for(counter = 0; static_cast<unsigned int>(counter) 
          < mcounters.size(); ++counter)
    {
      mcounters[counter] = counter*(newheight/mcounters.size());
      mstop[counter] = (counter+1)*(newheight/mcounters.size());
    }
    
    //set the last partition to the right size
    mstop[mcounters.size()-1] = newheight;

    projectPVFS(progress);


  }
  catch(...)
  {
    if (pmesh)
      delete pmesh;
    throw ProjectorException(PROJECTOR_ERROR_UNKOWN);

  }
}

//*****************************************************************
void PVFSProjector::projectPVFS(BaseProgress * progress)
  throw(ProjectorException)
{
  MPI_Status status;
  int msize(0), membersize(0), position(0); 
  long int buffersize(0);
  unsigned char * buffer(0);                     //the bufer for sending

  Stitcher * mystitch(0);                        //stitcher pointer
  long int beginofchunk(0), endofchunk(0);       //for passing to the slave
  long int chunkcounter(0);                      //for output
  long int ycounter(1);
  int chunkdif(maxchunk-minchunk);


  try
  {
                                                 //init the status progress
    if (progress)
    {
      std::strstream tempstream;
      
      tempstream << "Reprojecting " << newheight << " lines." << std::ends;
      tempstream.freeze(0);
      progress->init(tempstream.str(),
                     NULL,
                     "Done.",
                     newheight, 
                     29);
      progress->start();                         //start the progress
    }
  
    if (stitcher)                                //see if we want a stitcher  
    {
                                                 //creates the stitcher thread
      if (!(mystitch = new (std::nothrow) Stitcher(out)))
        
        throw std::bad_alloc();
    }
    
     //figure out the maximum buffer size based on the system;
    MPI_Pack_size(2, MPI_LONG, MPI_COMM_WORLD, &membersize);
    buffersize+=membersize;
    if(!slavelocal)
    {
      MPI_Pack_size(maxchunk*newwidth*spp, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD,
                    &membersize);
      buffersize+=membersize;
    }
    //ask for the buffer
    if (!(buffer = new unsigned char[buffersize]))
      throw std::bad_alloc();
    

                               
    if (sequencemethod == 2)                     //init the random number
      srand48(time(NULL));
    
    while (chunkcounter < newheight)
    {
      //do a blocking wait for any message.
      MPI_Recv(buffer, buffersize, MPI_PACKED, MPI_ANY_SOURCE,
               MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      
      //get the starting scanline
      beginofchunk = mcounters[membership[status.MPI_SOURCE]]; 
     
     
      //check termination
      if (beginofchunk < 0)
      {
        //terminate the slave
        chunkcounter += terminateSlave(status, 
                                       mystitch, buffer, buffersize);
      }
      else
      {
        //check the sequence method
        switch(sequencemethod)
        {
        case 0:
          endofchunk = beginofchunk + maxchunk-1;
          break;
        case 1:
          ++ycounter;
          if (ycounter >= sequencesize)
            ycounter = 0;
          endofchunk = beginofchunk + sequence[ycounter]-1;
          break;
        case 2:
          endofchunk = beginofchunk + static_cast<int>(drand48()*chunkdif
                                                       + minchunk) -1;
          break;
        }
        
        //check to see if this is the last chunk in the partition
        if (endofchunk >= mstop[membership[status.MPI_SOURCE]]-1)
        {
          endofchunk = mstop[membership[status.MPI_SOURCE]]-1;
          //reset the counter
          mcounters[membership[status.MPI_SOURCE]] = -1;
        }
        else
        {
          //update the counter
          mcounters[membership[status.MPI_SOURCE]]+=
            (endofchunk-beginofchunk)+1;
        }

        switch(status.MPI_TAG)
        {
        case SETUP_MSG:
          //pack the info in
          sendSlaveSetup(status.MPI_SOURCE);
         
           //add the first bit of work
          position = 0;
          MPI_Pack(&(beginofchunk), 1, MPI_LONG, buffer, buffersize,
                   &position, MPI_COMM_WORLD);
          MPI_Pack(&(endofchunk), 1, MPI_LONG, buffer, buffersize,
                   &position, MPI_COMM_WORLD);
          //send to the slave
          MPI_Send(buffer, position, MPI_PACKED, status.MPI_SOURCE,
                   WORK_MSG, MPI_COMM_WORLD);
          break;
        case WORK_MSG:
          
          MPI_Get_count(&status, MPI_PACKED, &msize);
          //unpack the scanline
          if(!slavelocal)
          {
            if (stitcher)
            {
              chunkcounter += sendStitcher(mystitch, buffer, msize);
            }
            else
              chunkcounter += unpackScanline(buffer, msize);
          }
          else
          {
            chunkcounter += unpackRawWork(buffer, msize); 
          }

            //pack the next work
          position = 0;
          MPI_Pack(&(beginofchunk), 1, MPI_LONG, buffer, buffersize,
                   &position, MPI_COMM_WORLD);
          MPI_Pack(&(endofchunk), 1, MPI_LONG, buffer, buffersize,
                   &position, MPI_COMM_WORLD);
          //send to the slave
          MPI_Send(buffer, position, MPI_PACKED, status.MPI_SOURCE,
                   WORK_MSG, MPI_COMM_WORLD);
          break;
        case ERROR_MSG:
        default:
          throw ProjectorException(PROJECTOR_ERROR_BADINPUT);
        }
      }
      
      //update the output
      if (progress && !((chunkcounter) % 11))
        progress->update(chunkcounter);

    }

    if (progress)
      progress->done();
    
    if (stitcher)
    {
      mystitch->wait();
      //remove the stitcher
      delete mystitch;
    }
    
    if(!slavelocal)
      writer.removeImage(0);                      //flush the output image

    out = NULL;
  }
  catch(...)
  {
    if (stitcher)
    {
      delete mystitch;                          //should stop the stitcher
      mystitch = NULL;
    }
    if(!slavelocal)
      writer.removeImage(0);
  }
  
}

//*************************************************************
long int PVFSProjector::terminateSlave(MPI_Status status,
                                       Stitcher * mystitch,  
                                       unsigned char * buffer,
                                       long int buffersize)
  throw(ProjectorException)
{
  int position(0);
  long int retvalue(0);                              //return written
  long int maxdif(0);                                //for membership repartion
  long int beginofchunk(0), endofchunk(0);           //chunksizes
  unsigned int counter(0);
  int msize(0);

  //check to see if the slave is has work to give
  if (status.MPI_TAG == WORK_MSG)
  {
    if(!slavelocal)
    {
      if (stitcher)
      {
        MPI_Get_count(&status, MPI_PACKED, &msize);
        retvalue = sendStitcher(mystitch, buffer, msize);
      }
      else
      {
        MPI_Get_count(&status, MPI_PACKED, &msize);
        retvalue = unpackScanline(buffer, msize);
      }
    }
    else
    {
      MPI_Get_count(&status, MPI_PACKED, &msize);
      retvalue = unpackRawWork(buffer, msize);

    } 
  }
  else if (status.MPI_TAG == ERROR_MSG)
  {
    throw ProjectorException(PROJECTOR_ERROR_BADINPUT);
    //should never happen :>
  }

  //see if we can change the slave nodes membership
  for(; counter < mcounters.size(); ++counter)
  {
    if (mcounters[counter] >= 0)
    {
      if (maxdif < (mstop[counter] - mcounters[counter]))
      {
	maxdif = (mstop[counter] - mcounters[counter]);
        membership[status.MPI_SOURCE] = counter; //change the membership
      }
    }
  }

  if (maxdif)
  {
    //found a good membership
    beginofchunk = mcounters[membership[status.MPI_SOURCE]];

    /**
     *TODO: Need to have a sequence checking function that
     *      will produce the next correct chunk no matter where
     *      it is called.
     **/
    endofchunk = beginofchunk + maxchunk-1;


    //check to see if this is the last chunk in the partition
    if (endofchunk >= mstop[membership[status.MPI_SOURCE]]-1)
    {
      endofchunk = mstop[membership[status.MPI_SOURCE]]-1;
      //reset the counter
      mcounters[membership[status.MPI_SOURCE]] = -1;
    }
    else
    {
      //update the counter
      mcounters[membership[status.MPI_SOURCE]]+=(endofchunk-beginofchunk)+1;
    }

    //pack the work and send the slave of to its new membership
    position = 0;
    MPI_Pack(&(beginofchunk), 1, MPI_LONG, buffer, buffersize,
             &position, MPI_COMM_WORLD);
    MPI_Pack(&(endofchunk), 1, MPI_LONG, buffer, buffersize,
             &position, MPI_COMM_WORLD);
          //send to the slave
    MPI_Send(buffer, position, MPI_PACKED, status.MPI_SOURCE,
             WORK_MSG, MPI_COMM_WORLD);
  }
  else
  {
    //the slave is done so it should get out of dodge
    MPI_Send(0, 0, MPI_PACKED, status.MPI_SOURCE,
              EXIT_MSG, MPI_COMM_WORLD);
  }

  return retvalue;

}

//*******************************************************************
void PVFSProjector::setUpRawOutput() throw(ProjectorException)
{
  double tp[6] = {0};
  double res[3] = {0};
  int ofiledesc(0);                         //output file descriptor
  pvfs_filestat metadata = {0, 0, 0, 0, 0}; //the meta data 
  try
  {
    //perform some checks
    if (!newwidth || !newheight || !infile)
      throw ProjectorException();
    
    //setup the meta data
    metadata.base = -1;   //start at default base
    metadata.pcount = -1; //go across as many nodes as possible?
    metadata.bsize = 0;   
    //set the stripe size to be a multiple of scanline chunk size
    metadata.ssize = maxchunk*newwidth*spp;
    
    
    //Create the output world file
    tp[3] = outRect.left;
    tp[4] = outRect.top;
    res[0] = newscale.x;
    res[1] = newscale.y;
    res[2] = 0;
    
    //write the world file to this directory for now
    if(!writer.writeTFW("out.TFW", tp, res))
    {
      throw ProjectionException();
    }
 
    writeImageMetrics("out.METRICS");
   
    //try to create the file
    if((ofiledesc = pvfs_open(outfile.c_str(), O_TRUNC|O_WRONLY|O_CREAT|O_META,
                              0777, &metadata, 0)) == -1)
    {
      throw ProjectorException();
    }
    
    
    //close the file 
    pvfs_close(ofiledesc);

    slavelocalpath = outfile; //give the slaves the output file name

  }
  catch(...)
  {
    throw ProjectorException(PROJECTOR_UNABLE_OUTPUT_SETUP);
  }
  
}

//****************************************************************
long int PVFSProjector::unpackRawWork(unsigned char * buffer, 
                                      long int buffersize) throw()
{
  long int scanlinenumber(0);
  long int endscanline(0); 
  int position(0);
  try
  {
    MPI_Unpack(buffer, buffersize, &position,
               &scanlinenumber, 1, MPI_LONG, MPI_COMM_WORLD);
    MPI_Unpack(buffer, buffersize, &position,
               &endscanline, 1, MPI_LONG, MPI_COMM_WORLD);

     return endscanline-scanlinenumber + 1;

  }
  catch(...)
  {
    return -1;
    //something really bad happened 
   
  }


}

//*****************************************************************
void PVFSProjector::writeImageMetrics(const std::string & filename) 
  throw(ProjectorException)
{
  std::ofstream out(filename.c_str());

  
  //output the image metrics
  out << "name: " << outfile << std::endl;
  out << "width: " << newwidth << std::endl;
  out << "height: " << newheight << std::endl;
  out << "spp: " << spp << std::endl;
  out << "bps: " << bps << std::endl;

  out.close();
}
  
  
