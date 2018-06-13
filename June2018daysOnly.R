 
  library(PBSmapping)     
  library(sas7bdat)
  library(date)
  library(crawl)
  library(circular)

  ####### I have found that R is too inefficient (as I program it) to ########
  ####### subset the data to first. last.observations each vr2-day so ########
  ####### some of these data are preprocessd in SAS see               ########
  ####### read_telemetry_update.sas last seen under SECR subd         ########
  
  ### root of where telemetry, bathymetry and coastal mapping data are located
 location <- 'C:/Users/MD56M/Dropbox/rdmarfin12'

 # rd_dat <- read.sas7bdat(paste(location,'MARFIN12/May2016/fnldata_donly_update.sas7bdat',sep='/')) 
  
  
 # library(RCurl)
  rd_dat <- read.sas7bdat('https://github.com/StPeteMdMurphy/RedDrum-Telemetry/raw/master/fnldata_donly_update.sas7bdat')
#  rd_dat <- read.sas7bdat(x) 
  
  
  
            #head(rd_dat)
  # define dates as date variables
    rd_dat[,'date'] <- mdy.date(rd_dat[,'month'],rd_dat[,'day'],rd_dat[,'year'])  # dim(rd_dat)
 
 
    #### matrix file of depth data
   egulf <- read.table('C:/Users/MD56M/Dropbox/rdmarfin12/bathy_egulf.txt',col.names=c('prex','y','prez'))           
        ## head(egulf)      
   ############### DETECTOR LOCATION DATA ########################
 lat <- tapply(rd_dat$Y_UTM,rd_dat$vr2,mean)   #dim(lon)
 lon <- tapply(rd_dat$X_UTM,rd_dat$vr2,mean)
   lon_ <-  t(t(lon));  vr2_no <- as.numeric(t(t(names(lon))));
   lat_ <-  t(t(lat));  
     #### survey event data in UTM's
    vr2_station <- cbind(vr2_no,lon_,lat_)
    colnames(vr2_station) <- c('EID','X','Y')
    sta <-as.EventData(vr2_station,projection='UTM',zone=17)  # station number and location in UTM (Universal Transverse Mercator)
     #### convert survey event data to lat/lon    
    sta_ <- convUL(sta,km=FALSE)
    vr2_lst <- sta_[,1]                      # final station list
    vr2_locs <- sta_[,(2:3)]                 # final station location
   ###############################################################
   
   ############## EVENT DATA Location in LL #####################################
    Event_locs <- cbind(1:length(rd_dat$Y_UTM),rd_dat$X_UTM,rd_dat$Y_UTM)
    colnames(Event_locs) <- c('EID','X','Y')
    elocs <-as.EventData(Event_locs,projection='UTM',zone=17)
     #### convert survey event data to lat/lon    
    elocs_ <- convUL(elocs,km=FALSE)
    rd_dat_ <- cbind(rd_dat,elocs_[,2],elocs_[,3])    # dim(rd_dat)     
     colnames(rd_dat_)[c(dim(rd_dat)[2]+1,dim(rd_dat)[2]+2)] <- c('GPS_W','GPS_N')           
    ################################################################
  #tag list                 
   tag_lst <- t(levels(as.factor(rd_dat_$tag)))
       
         #### VR2 symbol colors
       drt <- sta_[,1]/max(sta_[,1])
       vr2_col <- rgb(drt,seq(from=0,to=1,by=1/(length(sta_[,1])-1)),
                          seq(from=1,to=0,by=-1/(length(sta_[,1])-1)),1)      
                
        # what is the maximum number of different detection days possible, variable is temp
                       temp <- 0
                  for(u in 1:length(tag_lst))        # u=1                            
                          {
                 det_set <- rd_dat_[rd_dat_$tag==as.numeric(tag_lst[u]),]  #dim(rd_dat_[rd_dat_$tag=="13868",]) dim(rd_dat_)

                    if(dim(det_set)[1]>temp) { temp <- dim(det_set)[1]; longtag<-tag_lst[u]; }
                          }
   #####DEFINE VARIABLES NEEDED TOP DESCRIBE EACH FISH'S MOVEMENT #################################
   step_dist <- matrix(rep(NA,temp*length(tag_lst)),byrow=TRUE,nrow=length(tag_lst))
   step_brng <- matrix(rep(NA,temp*length(tag_lst)),byrow=TRUE,nrow=length(tag_lst))
   step_time <- matrix(rep(NA,temp*length(tag_lst)),byrow=TRUE,nrow=length(tag_lst))
   disp_dist <- matrix(rep(NA,temp*length(tag_lst)),byrow=TRUE,nrow=length(tag_lst))
   disp_brng <- matrix(rep(NA,temp*length(tag_lst)),byrow=TRUE,nrow=length(tag_lst))
   det_time  <- matrix(rep(NA,temp*length(tag_lst)),byrow=TRUE,nrow=length(tag_lst))
   det_vr2   <- matrix(rep(NA,temp*length(tag_lst)),byrow=TRUE,nrow=length(tag_lst))
   start_time <- rep(NA,temp)            
   class(det_time) <- "POSIXct"
   class(start_time) <- "POSIXct"
   # ###################################
   # for each fish tagged determine start distance jumps and bearings
        len_det_set <- rep(NA,length(tag_lst))
        rel_loc <- rep(NA,length(tag_lst))
        for(k in 1:length(tag_lst))   # k=1                                 
           {
           det_set <- rd_dat_[rd_dat_$tag==as.numeric(tag_lst[k]),]
           len_det_set[k] <- dim(det_set)[1]
              ##### need the detections to be in chronological order here ####
                det_set<-det_set[order(det_set$date),]
              
           start_time[k] <- det_set$date[1] 
           rel_loc[k] <- det_set$vr2[1]                       
                 if(dim(det_set)[1]>1)                        
                   {
          for(m in 1:(dim(det_set)[1]-1))     #m=2                          
            { 
        det_time[k,m] <- det_set$date[m] 
        det_vr2[k,m] <- det_set$vr2[m+1]    
       # distance and bearing from starting point
        # all to radians
         rad_lat <- det_set$GPS_N * pi/180
         rad_lon <- (det_set$GPS_W+360) * pi/180   # convert longitude to angle east of prime meridian
   ######## DISTANCE ########################################
      # from origin 
      del_lat <- rad_lat[1]-rad_lat[m+1]
      del_lon <- rad_lon[1]-rad_lon[m+1]
    intermed <- sin(del_lat/2)^2 + cos(rad_lat[m+1])*cos(rad_lat[1])*sin(del_lon/2)^2
      intermed2 <- 2*asin(min(1,sqrt(intermed)))
      disp_dist[k,m] <- 6371. * intermed2
      # this step
      del_lat <- rad_lat[m]-rad_lat[m+1]
      del_lon <- rad_lon[m]-rad_lon[m+1]
    intermed <- sin(del_lat/2)^2 + cos(rad_lat[m+1])*cos(rad_lat[m])*sin(del_lon/2)^2
      intermed2 <- 2*asin(min(1,sqrt(intermed)))
      step_dist[k,m] <- 6371. * intermed2
   ##########################################################
      
    ######## BEARING (output degrees) #########################
     # from origin
    if(disp_dist[k,m]>0) {
       temp_b <-  atan2( sin(rad_lon[m+1]-rad_lon[1])*cos(rad_lat[m+1]),
       cos(rad_lat[1])*sin(rad_lat[m+1])-sin(rad_lat[1])*cos(rad_lat[m+1])*cos(rad_lon[m+1]-rad_lon[1])  )*  
               180/pi
       disp_brng[k,m] <- (temp_b+360)%%360   # convert to true 0-360 degree bearing 
                       }
      # this step
    if(step_dist[k,m]>0) {      
       temp_b <-  atan2( sin(rad_lon[m+1]-rad_lon[m])*cos(rad_lat[m+1]),
       cos(rad_lat[m])*sin(rad_lat[m+1])-sin(rad_lat[m])*cos(rad_lat[m+1])*cos(rad_lon[m+1]-rad_lon[m])  )*  
               180/pi
        step_brng[k,m] <- (temp_b+360)%%360   # convert to true 0-360 degree bearing 
                       }     
   ############################################################
         }  # for loop within tag  m
                   } # if condition
                         } # for across tags k    
    ######################################################################
    ########### All data assembled ######################################
    ####################################################################                     
   ptm <- proc.time()   # start time the program run         
          ### how many observations for each fish?
           samplN <- matrix(rep(NA,3*length(tag_lst)),ncol=3,byrow=TRUE)
          for(k in 1:length(tag_lst))  
           {
           samplN[k,1] <- k 
           samplN[k,3] <- length(step_dist[k,!is.na(step_dist[k,])])
           samplN[k,2] <- tag_lst[,k] 
           }      
         #  t(t( samplN[samplN[,2]>10,]   ))
         ########################################
                
                                           
    ###### Run only for this one fish for now (sequence in tag number***
    ###########################                     
    fno <-39
    #################################  
 
    ##### getting working versions of distance and bearing for this fish without NA's
    ##### seems the data were not in temporal order so need to sort by time
     det_timeW <- det_time[fno,!is.na(det_time[fno,])] 
      det_timeW <- det_timeW[order(det_timeW)]    
                                 # length(det_timeW)
    step_distW<-step_dist[fno,!is.na(step_dist[fno,])]
      step_distW <- step_distW[order(det_timeW)]   
                               # length(step_distW)
     step_brngW<-step_brng[fno,1:length(step_distW)]
       step_brngW <- step_brngW[order(det_timeW)]  
                              #length(step_brngW)
     det_vr2W <- det_vr2[fno,!is.na(det_vr2[fno,])]
        det_vr2W <-  det_vr2W[order(det_timeW)]  
                                 # length(det_vr2W)
      vr2_colW <- vr2_col[match(det_vr2W,vr2_lst)] 
        vr2_colW <- vr2_colW[order(det_timeW)]  
    #######################################################
    #######################################################

   # this tracking starts at the first detection time and location
     startT <-  det_timeW[1]
     det_1 <- vr2_locs[det_vr2W[1]==vr2_lst,]
        rad_st_lat <- as.numeric(det_1[2])*pi/180
        rad_st_lon <- (360+as.numeric(det_1[1]))*pi/180 # convert longitude to angle east of prime meridian
  
    # what were observed locations within this array and hours after first detect?
         #   obs_det_ <- matrix( rep(NA,3*(length(det_timeW)-1)),byrow=FALSE,ncol=3 )
         #   obs_difft <- rep(NA,length(det_timeW))
            obs_dat=NULL
           for(i in 1:length(det_timeW))   #i=1   
             {
           tobs_det <- vr2_locs[det_vr2W[i]==vr2_lst,]
           tobs_difft <- round(difftime(det_timeW[i],det_timeW[1],units="days"),digits=0)
           obs_dat <- rbind(obs_dat,c(tobs_det$X,tobs_det$Y,tobs_difft))  # str(tobs_det)

        #   obs_det_[i,(1:2)] <- vr2_locs[det_vr2W[i]==vr2_lst,]
        #   obs_difft[i] <- round(difftime(det_timeW[i],det_timeW[1],units="hours"),digits=0)
             }
        #     obs_dat <- cbind(t(t(obs_difft))+1,obs_det_)

 ###### compress to only average location each day #######
   ave_lon <- t(t( tapply(obs_dat[,1],obs_dat[,3],mean)))
   ave_lat <- t(t( tapply(obs_dat[,2],obs_dat[,3],mean)))
   days <- t(t(as.numeric(rownames(ave_lat))))
    obs_dat <- cbind(ave_lon+360,ave_lat,days)   # convert longitude into angle east of prime meridian

    ##### total number of hours this fish was detected (time between first and last detection in hours) ... for walk steps #####
     difft <- round(difftime(det_timeW[length(det_timeW)],det_timeW[1],units="days"),digits=0)


    ############################################################################
    
  ##### correlated random walk model
  library(circular)
  
  # length of walk
  N <-difft     # this is an number of days detmination
  #includes first detection time and hourly intervals throughout observed time after that
  sim_times <- seq(as.POSIXlt(det_timeW[1]),as.POSIXlt(det_timeW[length(det_timeW)]),by="days")
  
  
  # color vector, symbol, and size for simulation symbols
   pt_col <- rep('black',N+1)
   pch_lst <- rep(20,N+1)
   pchsz_lst <- rep(0.2,N+1)
                                                     
    ########################################################
    #########################################################
    ########################################################
    #########################################################
    ########################################################
    #########################################################

  
  #### multiple iterations of this fishes simulated track
      numsims <- 1
     wbscale <- c(0.5)#,0.9,1.1,1.3) #1.0,2.0)       # weibel scalar
     cmu <- c(0)#,pi/2,pi,3*pi/2)  # direction (mu, radians)-- 0=east,pi/2=south,pi=west,3/2pi=north
     crho <- c(0.1)#,0.4,0.6)
     weib_shp <- c(2)
   
      # weib_shp <- c(3,2,1)  # gives skew probably 0.5 (high right skew) to 2 is best, at 5 slight left skew
      # wbscale <- c(0.5,0.8)       # weibel scalar -- scale of daily jump distances (0.5 is small -- 1.0 is medium-large)
      # cmu <- c(0* i,3/2*pi)  # direction (mu, radians) 0*pi=east,pi/2=South,pi=west,3/2*pi=north
      # crho <- c(0.0,0.01,0.02) # correlation of jumps...above 0.2 gets linear
      
      
      
           runs_ <- numsims*length(wbscale)*length(cmu)*length(crho)
      run_counter <- 0
    sum_distdiffs <- matrix(rep(NA,runs_*6),ncol=6,byrow=TRUE) 
    tst_temp <-1e09 
      for(d in 1:numsims)  #d=1
       {                                                
# make weibull distributed N number of step distances    N=1000
 par(mfrow=c(2,2))
      for(weibsc in c(wbscale))  
         {
  steps <- rweibull(N,weib_shp,weibsc);     # no obs., shape, and scale parameters   << Remember distance per time # N =500
  #     hist(rweibull(N,5,100))
  #     }
#     bearings
     for(directR in c(cmu))     #direction in radians #directR=cmu
       {
       for(r in c(crho))                              #r=crho
         {
  theta <- rwrappedcauchy(N,mu=directR,rho=r)+pi/2 
        # mu is direction originally E and counted counter-clockwise
        # adjusted to make clockwise with N=0, E=-pi/2, S=-pi, W=-3pi/2, etc.
  #rose.diag(theta,bins=24)

# cumulative angle (absolute orientation)
  Phi <- cumsum(theta)

# step distance components
  dX <- steps*cos(theta)   
  dY <- steps*sin(theta)           

# actual X-Y values
#  X<-cumsum(dX)                                                 
#  Y<-cumsum(dY)                                                 
                                                              
 # par(mfrow=c(2,1)) 
 # plot(X,Y) 
 # plot(lon,lat)
                                     
  
    # length(steps)  
         #### what are these actual locations? ###
         lat <- rep(NA,length(steps))
         lon <- rep(NA,length(steps))
         distST <-  rep(NA,length(steps))
         depth <-  rep(NA,length(steps))
     # from first det site
        distST[1] <- sqrt( dX[1]^2 + dY[1]^2 )   
     
    ##### new absolute location based on bearing and distance ####   
   lat[1] <- (180/pi)*( asin(sin(rad_st_lat)*cos(distST[1]/6371)+cos(rad_st_lat)*sin(distST[1]/6371)*cos(theta[1])) )
   lon[1] <- (180/pi)*( rad_st_lon+atan2(sin(theta[1])*sin(distST[1]/6371)*cos(rad_st_lat),
                                            cos(distST[1]/6371)-sin(rad_st_lat)*sin(lat[1])) )  
                                            
               hilon <- sort(egulf[lon[1]<egulf[,1],1],decreasing=FALSE)[1] #first longitude class east of location
               hilat<- sort(egulf[lat[1]<egulf[,2],2],decreasing=FALSE)[1] #first latitude north of location 
                depth[1] <- egulf[egulf[,1]==hilon & egulf[,2]==hilat,3]                         
   
         for(i in 2:length(steps))#i=2
          {
         distST[i] <- sqrt( dX[i]^2 + dY[i]^2 )       
            ### given these distances (distST) and bearings (Phi) ###
            ### what is location lat,lon?  Radius of earth = 6371 km
            
  lat[i] <- (180/pi)*(asin(sin(lat[i-1]*pi/180)*cos(distST[i]/6371)+cos(lat[i-1]*pi/180)*sin(distST[i]/6371)*cos(theta[i-1]))  )
  lon[i] <- (180/pi)*(lon[i-1]*pi/180+atan2(sin(theta[i-1])*sin(distST[i]/6371)*cos(lat[i-1]*pi/180),
                                 cos(distST[i]/6371)-sin(lat[i-1]*pi/180)*sin(lat[i]*pi/180)) )
   
   #how deep is this                              
   hilon <- sort(egulf[lon[i]<egulf[,1],1],decreasing=FALSE)[1] #first longitude class east of location
   hilat<- sort(egulf[lat[i]<egulf[,2],2],decreasing=FALSE)[1] #first latitude north of location 
   depth[i] <- egulf[egulf[,1]==hilon & egulf[,2]==hilat,3] 
   
     ## <<off the edge of known universe ##
    if (is.na(depth[i]))
     { 
      depth[i] <- depth[i-1] 
     lat[i] <- lat[i-1]
     lon[i] <- lon[i-1]
     }
     ######################################
   # too shallow or too deep set the next bearing to likely turn back.... 
       if (depth[i]< (-40))     #### too deep go east
         {
          theta[i+1] <- rwrappedcauchy(1,mu = 0,rho=0.99)+pi/2 
          }             ## rwrappedcauchy(1,mu = -pi/2,rho=0.99)+pi/2 
                          # rose.diag(rwrappedcauchy(N,mu=pi/2-pi,rho=0.8)+pi/2 ,bins=24) 
         if (depth[i]> (-5))     #### too shallow go west
         {
          theta[i+1] <- rwrappedcauchy(1,mu = -3*pi/2,rho=0.99)+pi/2 
          }
          
       # #### if on land do until the fish gets back in the water in one step #######
         iterd <- 0  
          while(depth[i]>0&iterd<=50)  # fish after step is on dry land
          {
          iterd<-iterd+1
          #choose a new bearing toward west
          newbrng <- rwrappedcauchy(1,mu=-3*pi/2,rho=0.9)+pi/2   
               # bearing degrees
               
                #  rose.diag(newbrng,bins=24)
          
  lat[i] <- (180/pi)*(asin(sin(lat[i-1]*pi/180)*cos(distST[i]/6371)+cos(lat[i-1]*pi/180)*sin(distST[i]/6371)*cos(newbrng))  )
  lon[i] <- (180/pi)*(lon[i-1]*pi/180+atan2(sin(newbrng)*sin(distST[i]/6371)*cos(lat[i-1]*pi/180),
                                 cos(distST[i]/6371)-sin(lat[i-1]*pi/180)*sin(lat[i]*pi/180)) )
     #how deep is this                              
   hilon <- sort(egulf[lon[i]<egulf[,1],1],decreasing=FALSE)[1] #first longitude class east of location
   hilat<- sort(egulf[lat[i]<egulf[,2],2],decreasing=FALSE)[1] #first latitude north of location 
   depth[i] <- egulf[egulf[,1]==hilon & egulf[,2]==hilat,3] 
           } 
           
       
            
        ############################################################################
                         # rose.diag(rwrappedcauchy(N,mu=-3*pi/2,rho=0.7)+pi/2 ,bins=24)       
          }  # end for loop determining new location
   
   ########### summary of simulation data ###############################     
   fnl_lat <- as.numeric(c(det_1$Y,lat))           
   fnl_lon <- as.numeric(c(det_1$X+360,lon))
   
     hilon <- sort(egulf[fnl_lon[1]<egulf[,1],1],decreasing=FALSE)[1] #first longitude class east of location
     hilat<- sort(egulf[fnl_lat[1]<egulf[,2],2],decreasing=FALSE)[1] #first latitude north of location 
                init_depth <- egulf[egulf[,1]==hilon & egulf[,2]==hilat,3]
   fnl_depth <- c(init_depth,depth)
   
        ### rounding in sim_times---add observation if so?
  if(length(sim_times)<length(fnl_lat)) sim_times <- c(sim_times,sim_times[length(sim_times)]+3600);
  if(length(theta)<length(sim_times)) theta <- c(theta,NA) 
  if(length(distST)<length(sim_times)) distST <- c(distST,NA) 
     
   sim_data <- cbind( t(t(sim_times)),t(t(fnl_lat)),t(t(fnl_lon)),
                      t(t(fnl_depth)),t(t(theta*180/pi)),t(t(distST)) )   #  length(distST)
                    
   #############################################################                     
   
  ##### what is the closest time in sim data to that in observed data? #####
    lookt <-  sort(sim_times,descending=FALSE) 
       #### look at lower and higher sim time step ##
          cntr <- rep(NA,length(det_timeW))
          cntr[1] <- 1
       for(k in 2:length(det_timeW))
         {
       laterT <- lookt[as.POSIXct(det_timeW[k])<lookt][1]
       earlierT <-  lookt[as.POSIXct(det_timeW[k])>lookt][1]
    ifelse(difftime(as.POSIXct(det_timeW[k]),as.POSIXct(earlierT))
             > difftime(as.POSIXct(laterT),as.POSIXct(det_timeW[k])),
       cntr[k] <- match(laterT,lookt), cntr[k] <- match(earlierT,lookt))  
  # for when the last observed detection is later than end of sim  
  if(is.na(cntr[k])) cntr[k] <- length(lookt)
          }
       
       ###### for the observation-simulation matches ####   
       pt_col[cntr]<- vr2_colW         # length(pt_col)
       pch_lst[cntr] <- 5 
       pchsz_lst[cntr] <- 0.5
   # HERE must determine distances between those common times in observed and simulation      
                                      distdiff <- rep(NA,dim(obs_dat)[1])
         for(h in 1:dim(obs_dat)[1])
           {
           #distance between common points
           rad_lat1 <- as.numeric(obs_dat[h,2])*pi/180
           rad_lon1 <- (as.numeric(obs_dat[h,1])+360)*pi/180
        #   rad_lat2 <- sim_data[cntr[h],2]*pi/180
        #   rad_lon2 <- sim_data[cntr[h],3]*pi/180
          rad_lat2 <- sim_data[(obs_dat[h,3]+1),2]*pi/180
          rad_lon2 <- sim_data[(obs_dat[h,3]+1),3]*pi/180
      del_lat <- rad_lat1-rad_lat2
      del_lon <- rad_lon1-rad_lon2
      intermed <- sin(del_lat/2)^2 + cos(rad_lat[m+1])*cos(rad_lat[m])*sin(del_lon/2)^2
      intermed2 <- 2*asin(min(1,sqrt(intermed)))
      distdiff[h] <- 6371. * intermed2
           }
     
            run_counter <- run_counter+1
          sum_distdiffs[run_counter,1] <- sum(distdiff)
          sum_distdiffs[run_counter,2] <- d # sim_number
          sum_distdiffs[run_counter,3] <- weibsc # weibel scale factor
          sum_distdiffs[run_counter,4] <- (directR+pi/2)*180/pi # compass heading
          sum_distdiffs[run_counter,5] <- r # correlation for compass heading 
          sum_distdiffs[run_counter,6] <- 1 # good run           
      if(iterd==50) sum_distdiffs[run_counter,6] <- 0;  # bad run   
                                         #min(sum_distdiffs[,1],na.rm=TRUE)
                                         #     hist(sum_distdiffs[sum_distdiffs[,2]==4&
                                         #                         sum_distdiffs[,3]==0.9,1])
             if (sum_distdiffs[run_counter,1] < tst_temp) 
                {
                  aside_sim <- sim_data
                  aside_obs <- obs_dat
                  tst_temp <- sum_distdiffs[run_counter,1]
                }
                
              } #correlation for random walk  
              } #####compass heading loop  
              
              } #### weibel scale parameter
                
                }  #### run simulations across d
     #              sum_distdiffs[sum_distdiffs[,1]==min(sum_distdiffs[,1]),]

 ##################################################################
 ##################################################################
 ##################################################################
 ##################################################################
 ##################################################################
 ##################################################################
 ##################################################################
 ##################################################################
    runtime <- proc.time()-ptm      #### how long did this take?
                                                               # warnings()
    par(mfrow=c(1,1))                                                                                       #hist(as.numeric(theta))              sin(theta[i])
# library(PBSmapping)      # read bathymetry file
       
#plotMap(egPoly,type="n")
clr <- .PBSclr()
 if(dev.cur() == 1) dev.new()
  #### telemetry event data
 walker <- data.frame(cbind(1:dim(aside_obs)[1],
                            as.numeric(aside_obs[,1]),
                            as.numeric(aside_obs[,2]),
                            as.numeric(aside_obs[,3])))

 colnames(walker) <- c('EID','X','Y','label')
  obs_ <-as.EventData(walker)
   attr(obs_,"projection") <-"LL"
  ### simulation data
   #### telemetry event data  head(aside_sim)
 stalker <- data.frame(cbind(t(t(1:dim(aside_sim)[1])),t(t(aside_sim[,3])),
                             t(t(aside_sim[,2]))),t(t(1:dim(aside_sim)[1])) )
 colnames(stalker) <- c('EID','X','Y','label')
  sim_ <-as.EventData(stalker)
   attr(sim_,"projection") <-"LL"

 sconnected <- as.PolySet(data.frame(PID=rep(1,dim(stalker)[1]),POS=1:dim(stalker)[1],X=stalker$X,Y=stalker$Y),projection="LL",zone=17)
 oconnected <- as.PolySet(data.frame(PID=rep(1,dim(walker)[1]),POS=1:dim(walker)[1],X=walker$X,Y=walker$Y),projection="LL",zone=17)

  ######################    head(sim_)    
                      
florida <- importGSHHS(paste(location,"gshhg-bin-2.2.3/gshhs_f.b",sep='/'),
  xlim=c(round(min(sim_[,2],obs_[,2])-0.1,1),round(max(sim_[,2],obs_[,2])+0.1,1)), 
  ylim=c(round(min(sim_[,3],obs_[,3])-0.1,1),round(max(sim_[,3],obs_[,3])+0.1,1)), n=15)
                             # longitude in degrees east of prime meridian
plotMap(florida, bg=clr$sea, col=clr$land)
#bathy
#addLines(egPoly,col=c('black','blue','blue','blue','blue'),lty=c(1,3,3,3,3))
   # addPoints(obs_,col=vr2_colW,pch=10)
       addPoints(obs_,pch=10,col='blue')
   
    addLines(oconnected,col='blue')

    addLines(sconnected,col='red')
    addPoints(sim_[1,],col="black",pch=21,cex=1.5)  # tag release spot   warnings()
 #   addPoints(sim_[2:dim(sim_)[1],],col=pt_col[2:length(pt_col)],
 #                                   pch=pch_lst[2:length(pch_lst)],
 #                                   cex=pchsz_lst[2:length(pchsz_lst)])
 #   addPoints(sim_[11:20,],col="orange",pch=7)
 #   addPoints(sim_[20:dim(sim_)[1],],col="green",pch=7)
#legend(x="bottomleft",bty="n",col=c('black','blue','blue','blue','blue'),
#         lwd=1,legend=as.character(isob),lty=c(1,3,3,3,3))
#title(for_title)
        
                      #  pch_lst <- rep(20,N+1)
                      #  pchsz_lst
    # head(sum_distdiffs)
    df <- data.frame(dev2 = sum_distdiffs[,1], rep = sum_distdiffs[,2],
                     wsf  = sum_distdiffs[,3], cmp = sum_distdiffs[,4],
                     cor  = sum_distdiffs[,5], flg = sum_distdiffs[,6])
   combos <-  unique(df[,c('wsf','cmp','cor')])
   
   # par(mfrow=c(4,4),mar=c(2,2,2,2))
   # 
   #     for(i in 1:dim(combos)[1]) # i=1
   #     {
   #       hist(df[df$wsf==combos[i,1]&df$cmp==combos[i,2]&df$cor==combos[i,3],1],breaks=8,
   #            main=paste(combos[i,1],combos[i,2],combos[i,3],sep='.'),xlab='')
   #     }
   # 
   #      testit <- lm(dev2~factor(wsf)+factor(cmp)+factor(cor),data=df)
   #      testit2 <- lm(dev2~factor(cmp)+factor(cor)+factor(wsf),data=df)
   #      anova(testit)
   #      anova(testit2)
   
    