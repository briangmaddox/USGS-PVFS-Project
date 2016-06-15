#ifndef MPIPROJECTORSLAVE_CPP
#define MPIPROJECTORSLAVE_CPP


#include "MpiProjectorSlave.h"
#include <unistd.h>

//********************************************************
MpiProjectorSlave::MpiProjectorSlave() : Projector(),
                                         slavelocal(false),
                                         mastertid(0), mytid(0),
                                         maxchunk(1)
{
}

//********************************************************
MpiProjectorSlave::~MpiProjectorSlave()
{
}

//********************************************************
bool MpiProjectorSlave::connect() throw()
{
  double x, y;                             //temp projecting vars
  int _x, _y;                              //actual image xy
  long int currenty, endy,xcounter,
    ycounter;                              //current line and counter
  unsigned char * scanline = NULL;         //output scanline
  unsigned char * buffer = NULL;           //the buffer to send back
  const unsigned char * inscanline = NULL; //input scaline
  int sppcounter;                          // spp counter
  PmeshLib::ProjectionMesh * pmesh = NULL; //projection mesh
  double xscaleinv, yscaleinv;
  unsigned char * sendb(0);                //the send buffer
  int sendbsize(0);                        //the send buffer size
  MPI_Status status;                       //mpi status
  int msize(0),                            //the message size
    position;                              //for MPI unpacking
  
  try
  {
   
    //get my rank
    MPI_Comm_rank(MPI_COMM_WORLD, &mytid);
    
    //hope mpi lets you send zero length messages
    MPI_Send(0, 0, MPI_PACKED, 0,
                 SETUP_MSG, MPI_COMM_WORLD);
  
    unpackSetup();                      //unpack the setup info

    if (slavelocal)                     //check the slavelocal
    {
      return storelocal();            
    }

  
    if (pmesh)
    {
      delete pmesh;
      pmesh = NULL;
    }
     
    pmesh = setupReversePmesh();        //setup the reverse pmesh
     
    
    //create the buffer to be at least as big as the 
    //maximum chunksize
    if (!(buffer = new (std::nothrow) unsigned char 
          [(maxchunk)*newwidth*spp]))
      throw std::bad_alloc();
    

    //get the send back size
     //figure out the maximum buffer size based on the system;
    MPI_Pack_size(2, MPI_LONG, MPI_COMM_WORLD, &msize);
    sendbsize += msize;
    MPI_Pack_size(maxchunk*newwidth*spp, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD,
                  &msize);
    sendbsize+= msize;
  
    //create the mpi buffer
    if (!(sendb = new unsigned char [sendbsize]))
      throw std::bad_alloc();

   

    //get the first bit of work
    MPI_Recv(sendb, sendbsize, MPI_PACKED, MPI_ANY_SOURCE,
               MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    
    //calcuate the inversers to do multiplaction instead of division
    xscaleinv = 1.0/oldscale.x;
    yscaleinv = 1.0/oldscale.y;


    //proccess messages
    while (status.MPI_TAG != EXIT_MSG)
    {
      //get the size of the message
      MPI_Get_count(&status, MPI_PACKED, &msize);

      position = 0;
      MPI_Unpack(sendb, msize, &position, &currenty, 1, MPI_LONG,
                 MPI_COMM_WORLD);
      MPI_Unpack(sendb, msize, &position, &endy, 1, MPI_LONG,
                 MPI_COMM_WORLD);
      
      //now just build the return message
      position = 0;
      MPI_Pack(&currenty, 1, MPI_LONG, sendb, sendbsize, &position,
               MPI_COMM_WORLD);
      MPI_Pack(&endy, 1, MPI_LONG, sendb, sendbsize, &position,
               MPI_COMM_WORLD);

             
      for (ycounter = currenty; ycounter <= endy; ++ycounter)
      {
        scanline = &(buffer[newwidth*spp*(ycounter-currenty)]);
        
        //reproject the line
        for (xcounter = 0; xcounter < newwidth; ++xcounter)
        {
          x = outRect.left + newscale.x * xcounter;
          y = outRect.top  - newscale.y * ycounter;
          
          //now get the 
          //reverse projected value
          if (pmesh)
	  {
	    pmesh->projectPoint(x, y);
	  }        
          else
          {
            toprojection->projectToGeo(x, y, y, x);
            fromprojection->projectFromGeo(y, x, x, y);
          }
          _x = static_cast<long int>((x - inRect.left)* 
                                     (xscaleinv) + 0.5);
          _y = static_cast<long int>((inRect.top - y) * 
                                     (yscaleinv) + 0.5);
          
          if ((_x >= oldwidth) || (_x < 0) || (_y >= oldheight) 
              || (_y < 0))
          {
            for (sppcounter = 0; sppcounter < spp; sppcounter++)
              scanline[xcounter*spp + sppcounter] = 0;
          }
          else
          {
            inscanline = cache->getRawScanline(_y);
            for (sppcounter = 0; sppcounter < spp; sppcounter++)
              scanline[xcounter*spp + sppcounter ] = 
                inscanline[_x*spp + sppcounter];
          }
       
         
	}

      }
     
      //pack this chunk into the buffer
      MPI_Pack(buffer, (endy-currenty + 1)*newwidth*spp, 
               MPI_UNSIGNED_CHAR, sendb, sendbsize, &position,
               MPI_COMM_WORLD);
      
      

      //send the entire chunk back the the master
      MPI_Send(sendb, position, MPI_PACKED, 0,
                 WORK_MSG, MPI_COMM_WORLD);
      
      //reset the scanline
      scanline = NULL;

      //get the next work
      MPI_Recv(sendb, sendbsize, MPI_PACKED, MPI_ANY_SOURCE,
             MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    
      
    }
    
    delete [] buffer;
    scanline = NULL;
    delete pmesh;
    delete toprojection;
    toprojection = NULL;
    pmesh = NULL;
    return true;
  }
  catch(...)
  {
     //set a error to the master
    MPI_Send(0, 0, MPI_PACKED, 0,
             ERROR_MSG, MPI_COMM_WORLD);
    delete pmesh;
    delete toprojection;
    toprojection = NULL;
    pmesh = NULL;
    delete scanline;
    return false;
  }
}

//*********************************************************
bool MpiProjectorSlave::storelocal() throw()
{
  double x, y;                             //temp projecting vars
  int _x, _y;                              //actual image xy
  long int currenty, endy,xcounter,
    ycounter;                              //current line and counter
  unsigned char * scanline = NULL;         //output scanline
  unsigned char * buffer = NULL;           //the buffer to send back
  const unsigned char * inscanline = NULL; //input scaline
  int sppcounter;                          // spp counter
  PmeshLib::ProjectionMesh * pmesh = NULL; //projection mesh
  double xscaleinv, yscaleinv;
  unsigned char * sendb(0);                //the send buffer
  int sendbsize(0);                        //the send buffer size
  MPI_Status status;                       //mpi status
  int msize(0),                            //the message size
    position;                              //for MPI unpacking
  int ofiledesc(0);                        //output file desc
  

  try
  {
   
    if((ofiledesc = pvfs_open(basepath.c_str(), O_WRONLY, 0777, 0, 0)) == -1)
    {
      //throw out
      throw std::bad_alloc();
    }
    
    

  
    if (pmesh)
    {
      delete pmesh;
      pmesh = NULL;
    }
     
    pmesh = setupReversePmesh();        //setup the reverse pmesh
     
    
    //create the buffer to be at least as big as the 
    //maximum chunksize
    if (!(buffer = new (std::nothrow) unsigned char 
          [(maxchunk)*newwidth*spp]))
      throw std::bad_alloc();
    

    //get the send back size
     //figure out the maximum buffer size based on the system;
    MPI_Pack_size(2, MPI_LONG, MPI_COMM_WORLD, &msize);
    sendbsize += msize;

    /*MPI_Pack_size(maxchunk*newwidth*spp, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD,
                  &msize);
    sendbsize+= msize;
    */
    
    //create the mpi buffer
    if (!(sendb = new unsigned char [sendbsize]))
      throw std::bad_alloc();

    //get the first bit of work
    MPI_Recv(sendb, sendbsize, MPI_PACKED, MPI_ANY_SOURCE,
               MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
    //calcuate the inversers to do multiplaction instead of division
    xscaleinv = 1.0/oldscale.x;
    yscaleinv = 1.0/oldscale.y;


    //proccess messages
    while (status.MPI_TAG != EXIT_MSG)
    {
      //get the size of the message
      MPI_Get_count(&status, MPI_PACKED, &msize);

      position = 0;
      MPI_Unpack(sendb, msize, &position, &currenty, 1, MPI_LONG,
                 MPI_COMM_WORLD);
      MPI_Unpack(sendb, msize, &position, &endy, 1, MPI_LONG,
                 MPI_COMM_WORLD);
      
      //now just build the return message
      position = 0;
      MPI_Pack(&currenty, 1, MPI_LONG, sendb, sendbsize, &position,
               MPI_COMM_WORLD);
      MPI_Pack(&endy, 1, MPI_LONG, sendb, sendbsize, &position,
               MPI_COMM_WORLD);

             
      for (ycounter = currenty; ycounter <= endy; ++ycounter)
      {
        scanline = &(buffer[newwidth*spp*(ycounter-currenty)]);
        
        //reproject the line
        for (xcounter = 0; xcounter < newwidth; ++xcounter)
        {
          x = outRect.left + newscale.x * xcounter;
          y = outRect.top  - newscale.y * ycounter;
          
          //now get the 
          //reverse projected value
          if (pmesh)
	  {
	    pmesh->projectPoint(x, y);
	  }        
          else
          {
            toprojection->projectToGeo(x, y, y, x);
            fromprojection->projectFromGeo(y, x, x, y);
          }
          _x = static_cast<long int>((x - inRect.left)* 
                                     (xscaleinv) + 0.5);
          _y = static_cast<long int>((inRect.top - y) * 
                                     (yscaleinv) + 0.5);
          
          if ((_x >= oldwidth) || (_x < 0) || (_y >= oldheight) 
              || (_y < 0))
          {
            for (sppcounter = 0; sppcounter < spp; sppcounter++)
              scanline[xcounter*spp + sppcounter] = 0;
          }
          else
          {
            inscanline = cache->getRawScanline(_y);
            for (sppcounter = 0; sppcounter < spp; sppcounter++)
              scanline[xcounter*spp + sppcounter ] = 
                inscanline[_x*spp + sppcounter];
          }
       
         
	}

      }
     
      //seek to the right position in the file....
      //Maybe this should be a lseek64??
      if(pvfs_lseek(ofiledesc, currenty*newwidth*spp, SEEK_SET) == -1)
      {
        //throw out of it
        throw std::bad_alloc();
      }
     
      //write the buffer to pvfs
      if(pvfs_write(ofiledesc, reinterpret_cast<char*>(buffer), 
                    (endy-currenty + 1)*newwidth *spp) !=
         (endy-currenty + 1)*newwidth*spp)
      {
        //throw out
        throw std::bad_alloc();
      }
 
      //send the entire chunk back the the master
      MPI_Send(sendb, position, MPI_PACKED, 0,
                 WORK_MSG, MPI_COMM_WORLD);
      
      //reset the scanline
      scanline = NULL;

      //get the next work
      MPI_Recv(sendb, sendbsize, MPI_PACKED, MPI_ANY_SOURCE,
             MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    
      
    }

    //close the pvfs file
    pvfs_close(ofiledesc);
    
    delete [] buffer;
    scanline = NULL;
    delete pmesh;
    delete toprojection;
    toprojection = NULL;
    pmesh = NULL;
    return true;
  }
  catch(...)
  {
    pvfs_close(ofiledesc);
    //set a error to the master
    MPI_Send(0, 0, MPI_PACKED, 0,
             ERROR_MSG, MPI_COMM_WORLD);
    delete pmesh;
    delete toprojection;
    toprojection = NULL;
    pmesh = NULL;
    delete scanline;
    return false;
  }

 
}
  



//*********************************************************
void MpiProjectorSlave::unpackSetup() throw()
{
 
  char tempbuffer[100];
  int temp;
  unsigned char * buf(0);
  int bufsize(0), tempsize(0);
  int position(0);
  MPI_Status status;
  std::string inputfilename;

  try
  {
    //calculate the buffersize
    MPI_Pack_size(200, MPI_CHAR, MPI_COMM_WORLD, &tempsize);
    bufsize += tempsize;
    MPI_Pack_size(2, MPI_LONG, MPI_COMM_WORLD, &tempsize);
    bufsize += tempsize;
    MPI_Pack_size(23, MPI_DOUBLE, MPI_COMM_WORLD, &tempsize);
    bufsize += tempsize;
    MPI_Pack_size(8, MPI_INT, MPI_COMM_WORLD, &tempsize);
    bufsize += tempsize;
    
    //create the buffer
    if (!(buf = new (std::nothrow) unsigned char[bufsize]))
      throw std::bad_alloc();
    
    //receive the setup
    MPI_Recv(buf, bufsize, MPI_PACKED, MPI_ANY_SOURCE,
               MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    position = 0;

    //use the actual size
    MPI_Get_count(&status, MPI_PACKED, &bufsize);
      
    //unpack the filename
    MPI_Unpack(buf, bufsize, &position, tempbuffer, 100, MPI_CHAR,
               MPI_COMM_WORLD);
    inputfilename = tempbuffer;


    //image metrics
    MPI_Unpack(buf, bufsize, &position, 
               &newheight, 1, MPI_LONG, MPI_COMM_WORLD);
    MPI_Unpack(buf, bufsize, &position, &newwidth, 1, MPI_LONG,
            MPI_COMM_WORLD);
    MPI_Unpack(buf, bufsize, &position, &newscale.x, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
    MPI_Unpack(buf, bufsize, &position, &newscale.y, 1, MPI_DOUBLE,
            MPI_COMM_WORLD);
    MPI_Unpack(buf, bufsize, &position, &outRect.left, 1, MPI_DOUBLE,
            MPI_COMM_WORLD);
    MPI_Unpack(buf, bufsize, &position, &outRect.top, 1, MPI_DOUBLE,
            MPI_COMM_WORLD);
    MPI_Unpack(buf, bufsize, &position, &outRect.bottom, 1, MPI_DOUBLE,
            MPI_COMM_WORLD);
    MPI_Unpack(buf, bufsize, &position, &outRect.right, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  
    //unpack slave local
    MPI_Unpack(buf, bufsize, &position, &temp, 1, MPI_INT,
              MPI_COMM_WORLD);
    slavelocal = temp;
   
    //unpack the max chunksize
    MPI_Unpack(buf, bufsize, &position, &maxchunk, 1, MPI_INT,
            MPI_COMM_WORLD);

    MPI_Unpack(buf, bufsize, &position, tempbuffer, 100, MPI_CHAR,
             MPI_COMM_WORLD);
    
    basepath = tempbuffer;
    
    //unpack pmesh info
    MPI_Unpack(buf, bufsize, &position, &pmeshsize, 1, MPI_INT,
            MPI_COMM_WORLD);
    MPI_Unpack(buf, bufsize, &position, &pmeshname, 1, MPI_INT,
           MPI_COMM_WORLD);
    
    //pack the projection parameters
    MPI_Unpack(buf, bufsize, &position, 
               reinterpret_cast<int *>(&Params.projtype), 1, MPI_INT,
             MPI_COMM_WORLD);
    
    MPI_Unpack(buf, bufsize, &position, 
               reinterpret_cast<int *>(&Params.datum), 1, MPI_INT,
             MPI_COMM_WORLD);

    MPI_Unpack(buf, bufsize, &position, 
               reinterpret_cast<int *>(&Params.unit), 1, MPI_INT,
             MPI_COMM_WORLD);
   
    MPI_Unpack(buf, bufsize, &position, &Params.StdParallel1, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
    
    MPI_Unpack(buf, bufsize, &position, &Params.StdParallel2, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);

    MPI_Unpack(buf, bufsize, &position, &Params.NatOriginLong, 1, MPI_DOUBLE,
              MPI_COMM_WORLD);

    MPI_Unpack(buf, bufsize, &position, &Params.NatOriginLat, 1, MPI_DOUBLE,
              MPI_COMM_WORLD);
    
    MPI_Unpack(buf, bufsize, &position, &Params.FalseOriginLong, 1, MPI_DOUBLE,
              MPI_COMM_WORLD);
    
    MPI_Unpack(buf, bufsize, &position, &Params.FalseOriginLat, 1, MPI_DOUBLE,
              MPI_COMM_WORLD);
    
    MPI_Unpack(buf, bufsize, &position, &Params.FalseOriginEasting, 1, 
               MPI_DOUBLE,
             MPI_COMM_WORLD);   
    
    MPI_Unpack(buf, bufsize, &position, &Params.FalseOriginNorthing, 1, 
               MPI_DOUBLE,
             MPI_COMM_WORLD);
    
    MPI_Unpack(buf, bufsize, &position, &Params.CenterLong, 1, 
               MPI_DOUBLE,
              MPI_COMM_WORLD);
    
    MPI_Unpack(buf, bufsize, &position, &Params.CenterLat, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
    
    MPI_Unpack(buf, bufsize, &position, &Params.CenterEasting, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
    
    MPI_Unpack(buf, bufsize, &position, &Params.CenterNorthing, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);

    MPI_Unpack(buf, bufsize, &position, &Params.ScaleAtNatOrigin, 1, 
               MPI_DOUBLE,
              MPI_COMM_WORLD);
    
    MPI_Unpack(buf, bufsize, &position, &Params.AzimuthAngle, 1, MPI_DOUBLE,
              MPI_COMM_WORLD);
    
    MPI_Unpack(buf, bufsize, &position, &Params.StraightVertPoleLong, 1, 
               MPI_DOUBLE,
              MPI_COMM_WORLD);
    
    MPI_Unpack(buf, bufsize, &position, &Params.zone, 1, MPI_INT,
             MPI_COMM_WORLD);
    
    MPI_Unpack(buf, bufsize, &position, &Params.FalseEasting, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
    
    MPI_Unpack(buf, bufsize, &position, &Params.FalseNorthing, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
      
    Projector::setInputFile(inputfilename);//setup cache, input image metrics, 
                                           //input projection
    
    toprojection = SetProjection(Params); //get the to projection
    
    delete [] buf; //done with buffer
    
  }
  catch(...)
  {
    //set a error to the master
    MPI_Send(0, 0, MPI_PACKED, 0,
                 ERROR_MSG, MPI_COMM_WORLD);
    delete [] buf;
    delete toprojection;
    toprojection = NULL;
  }
}

#endif




