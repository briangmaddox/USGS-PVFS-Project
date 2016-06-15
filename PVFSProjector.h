/**
 *This is a extension onto the pvm project class
 *that is designed to take advantage of pvfs partitioning
 *By Chris Bilderback
 **/

#ifndef PVFSPROJECTOR_H
#define PVFSPROJECTOR_H

//pvfs includes for pvfs writing
#include <pvfs.h>
#include <pvfs_proto.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <malloc.h>

#include "MpiProjector.h"
#include <vector>
#include <hash_map>
#include <fstream>




//PVFSProject class partitions slaves to
//work on only a specific portion of the file
//which can reduce network traffic as less slaves
//are talking to the master to get their data.
//Class also has the ablity to write to a raw file
//on a pvfs partition.
class PVFSProjector : public MpiProjector
{
 public:

  //Constructor and destructor
  PVFSProjector();
  ~PVFSProjector();

  //This function indicates the number of pvfs partitions that
  //is being worked with. (Default is 2)
  void setPartitionNumber(const unsigned int & inPart) throw();
  unsigned int getPartitionNumber() const throw();

  //This function sets whether the pvfs project will write a 
  //raw file output file to a pvfs partition or not.
  //Default value (when projector is constructed) is false.)
  void writeToPVFS(bool inWriteRaw) throw();

  //main function which runs the projection
  virtual void 
    project(BaseProgress * progress = NULL) 
    throw(ProjectorException);

 protected:
  
  //This function actually projects the file
  void projectPVFS(BaseProgress * progress = NULL)
    throw(ProjectorException);

  //This function tells the slave to terminate
  long int PVFSProjector::terminateSlave(MPI_Status status,
                          Stitcher * mystitch,  
                          unsigned char * buffer,
                          long int buffersize) throw(ProjectorException);

  //Function that sets up the raw output file (if desired)
  void setUpRawOutput() throw(ProjectorException);

  //Function to unpack raw output from the slave
  long int unpackRawWork(unsigned char * buffer, 
                         long int buffersize) throw();

  //Function to write a image metrics file (width, height, spp, bpp)
  //to disk (for raw file)
  void writeImageMetrics(const std::string & filename) 
    throw(ProjectorException);


  std::vector<long int> mcounters;    //membership counters
  std::vector<long int> mstop;        //where to stop each membership
  std::hash_map<int, unsigned int> 
    membership;                       //nodal membership map.
};

#endif
