# Movement Trajectory Simulations from Abrahms et al. 2016

########################################################################
########################### Central Place Forager ######################

# Set starting coordinates, can be offset if spatially explicit
x_init <- 0
y_init <- 0
n.steps <- 3600
HR <- 10000
core <- 1000
dist <- 0
ang.temp <- 0

# Latitude: 111,111 meters = 1 degree
# Longitude: 111,111 meters * cos(latitude) = 1 degree

df.central <- data.frame(matrix(n.steps, 4))

for (i in 1:n.steps){
  
  # Function normally if distance from origin is less than home range size 
  if (dist < HR && dist < core) {
    
    step <- rgamma(1, shape = 2, rate = 2) * 500
    
    # Uniform distribution for direction
    heading <- runif(1,0.0001,360)
    
    # Determine the angle to perform trig
    angle <- 0 
    
    if (heading > 0 && heading <= 90) {
      angle <- (90 - heading)
      #print(angle)
    } else if (heading > 90 && heading <= 180) {
      angle <- (heading - 90)
      #print(angle)
    } else if (heading > 180 && heading <= 270) {
      angle <- (270 - heading)
      #print(angle)
    } else {
      angle <- (heading - 270)
      #print(angle)
    }
    
    # Use sine to determine the movement in y (latitude)
    rad_y <- angle*0.0174532925
    y_change <- sin(rad_y)*step
    y_change
    
    # Use cosine to determine the movement in x (longitude)
    rad_x <- angle*0.0174532925
    x_change <- cos(rad_x)*step
    x_change
    
    #Determine whether to add or subtract new value based on heading
    
    if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
      y_new <- y_init + y_change
      #print(y_new)
    } else {
      y_new <- y_init - y_change
      #print(y_new)
    }
    
    if (heading > 0 && heading <= 180) {
      x_new <- x_init + x_change
      #print(x_new)
    } else {
      x_new <- x_init - x_change
      #print(x_new)
    }
    
    y_init <- y_new
    x_init <- x_new
    
    dist <- sqrt(x_new^2 + y_new^2)
    ang.temp <- atan2(y_new, x_new)
    ang.temp <- ang.temp/0.0174532925
    if (ang.temp > 0 && ang.temp <= 90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else if (ang.temp > 90 && ang.temp <= 180) {
      ang.temp <- (360 - ang.temp)
      #print(angle)
    } else if (ang.temp < 0 && ang.temp >= -90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    }
    
    df.central[i,1] <- x_init #lat
    df.central[i,2] <- y_init #long
    df.central[i,3] <- dist
    df.central[i,4] <- ang.temp
    
  }
  
  # If distance is greater than home range, turn around from relative position  
  else if (dist > HR && dist > core) { 
    
    if (ang.temp < 180){
      ang.return <- 180 + ang.temp
    } else {
      ang.return <- ang.temp - 180
    }
    
    # As maintain heading towards origin if outside of core area
    
    step <- rgamma(1, shape = 2, rate = 2) * 1000
    
    heading <- rnorm(1,ang.return,45)
    
    angle <- 0 
    
    if (heading > 0 && heading <= 90) {
      angle <- (90 - heading)
      #print(angle)
    } else if (heading > 90 && heading <= 180) {
      angle <- (heading - 90)
      #print(angle)
    } else if (heading > 180 && heading <= 270) {
      angle <- (270 - heading)
      #print(angle)
    } else {
      angle <- (heading - 270)
      #print(angle)
    }
    
    # Use sine to determine the movement in y (latitude)
    rad_y <- angle*0.0174532925
    y_change <- sin(rad_y)*step
    y_change
    
    # Use cosine to determine the movement in x (longitude)
    rad_x <- angle*0.0174532925
    x_change <- cos(rad_x)*step
    x_change
    
    #Determine whether to add or subtract new value based on heading
    
    if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
      y_new <- y_init + y_change
      #print(y_new)
    } else {
      y_new <- y_init - y_change
      #print(y_new)
    }
    
    if (heading > 0 && heading <= 180) {
      x_new <- x_init + x_change
      #print(x_new)
    } else {
      x_new <- x_init - x_change
      #print(x_new)
    }
    
    y_init <- y_new
    x_init <- x_new
    
    dist <- sqrt(x_new^2 + y_new^2)
    ang.temp <- atan2(y_new, x_new)
    ang.temp <- ang.temp/0.0174532925
    if (ang.temp > 0 && ang.temp <= 90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else if (ang.temp > 90 && ang.temp <= 180) {
      ang.temp <- (360 - ang.temp)
      #print(angle)
    } else if (ang.temp < 0 && ang.temp >= -90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    }
    
    df.central[i,1] <- x_init #lat
    df.central[i,2] <- y_init #long
    df.central[i,3] <- dist
    df.central[i,4] <- ang.temp
  }
  
  # If inside of the core area, return to normal functioning
  else if (dist < HR && dist > core) {if (ang.temp < 180){
    ang.return <- 180 + ang.temp
  } else {
    ang.return <- ang.temp - 180
  }
    
    # As maintain heading towards origin if outside of core area
    
    step <- rgamma(1, shape = 2, rate = 2) * 750
    
    heading <- rnorm(1,ang.return,45)
    
    angle <- 0 
    
    if (heading > 0 && heading <= 90) {
      angle <- (90 - heading)
      #print(angle)
    } else if (heading > 90 && heading <= 180) {
      angle <- (heading - 90)
      #print(angle)
    } else if (heading > 180 && heading <= 270) {
      angle <- (270 - heading)
      #print(angle)
    } else {
      angle <- (heading - 270)
      #print(angle)
    }
    
    # Use sine to determine the movement in y (latitude)
    rad_y <- angle*0.0174532925
    y_change <- sin(rad_y)*step
    y_change
    
    # Use cosine to determine the movement in x (longitude)
    rad_x <- angle*0.0174532925
    x_change <- cos(rad_x)*step
    x_change
    
    #Determine whether to add or subtract new value based on heading
    
    if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
      y_new <- y_init + y_change
      #print(y_new)
    } else {
      y_new <- y_init - y_change
      #print(y_new)
    }
    
    if (heading > 0 && heading <= 180) {
      x_new <- x_init + x_change
      #print(x_new)
    } else {
      x_new <- x_init - x_change
      #print(x_new)
    }
    
    y_init <- y_new
    x_init <- x_new
    
    dist <- sqrt(x_new^2 + y_new^2)
    ang.temp <- atan2(y_new, x_new)
    ang.temp <- ang.temp/0.0174532925
    if (ang.temp > 0 && ang.temp <= 90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else if (ang.temp > 90 && ang.temp <= 180) {
      ang.temp <- (360 - ang.temp)
      #print(angle)
    } else if (ang.temp < 0 && ang.temp >= -90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    }
    
    df.central[i,1] <- x_init #lat
    df.central[i,2] <- y_init #long
    df.central[i,3] <- dist
    df.central[i,4] <- ang.temp
  }
}

start.date <- '2015-01-01' 
start.time <- '00:00:00' 
interval <- 60 

increment.mins <- interval * 60 
x <- paste(start.date, start.time) 

for(i in 1:n.steps) { 
  df.central[i,5] <- as.POSIXct(strptime(x, "%Y-%m-%d %H:%M:%S")  + i*increment.mins)
} 

class(df.central[,5]) <- 'POSIXct'

plot(df.central[,1], df.central[,2], type='b')
hist(df.central[,3])

colnames(df.central) <- c("Longitude","Latitude","Dist_from_Center","Angle_from Center","Datetime")

########################################################################
################################## Nomads ##############################
# Set starting coordinates, can be offset if spatially explicit
x_init <- 0
y_init <- 0
n.steps <- 3600
dist <- 0
ang.temp <- 0

state.df <- data.frame(matrix(0,n.steps,2))
df.nomad <- data.frame(matrix(0,n.steps,4))

for (i in 1:n.steps) {
  
  if (i > 1) {
    prev.state <- state.df[i-1,2]
    rand <- runif(1,0,1)
    
    if (rand >= 0.05 && prev.state == 'Forage' ) {
      state.df[i,1] <- rand
      state.df[i,2] <- 'Forage'
    } else if (rand < 0.05 && prev.state == 'Forage') {
      state.df[i,1] <- rand
      state.df[i,2] <- 'Search'
    } else if (rand >= 0.05 && prev.state == 'Search') {
      state.df[i,1] <- rand
      state.df[i,2] <- 'Search'
    } else if (rand < 0.05 && prev.state == 'Search') {
      state.df[i,1] <- rand
      state.df[i,2] <- 'Forage'
    }
    
  } else {
    rand <- runif(1,0,1)
    if (rand >= 1.0) {
      state.df[i,1] <- rand
      state.df[i,2] <- 'Forage'
    } else {
      state.df[i,1] <- rand
      state.df[i,2] <- 'Search'
    }
  }
}

for (i in 1:n.steps) {
  if (i > 1) {
    if (state.df[i,2] == 'Forage') {
      
      step <- rgamma(1, shape = 2, rate = 2) * 100
      
      # Uniform distribution for direction
      heading <- runif(1,0.0001,360)
      
      # Determine the angle to perform trig
      angle <- 0 
      
      if (heading > 0 && heading <= 90) {
        angle <- (90 - heading)
        #print(angle)
      } else if (heading > 90 && heading <= 180) {
        angle <- (heading - 90)
        #print(angle)
      } else if (heading > 180 && heading <= 270) {
        angle <- (270 - heading)
        #print(angle)
      } else {
        angle <- (heading - 270)
        #print(angle)
      }
      
      # Use sine to determine the movement in y (latitude)
      rad_y <- angle*0.0174532925
      y_change <- sin(rad_y)*step
      y_change
      
      # Use cosine to determine the movement in x (longitude)
      rad_x <- angle*0.0174532925
      x_change <- cos(rad_x)*step
      x_change
      
      #Determine whether to add or subtract new value based on heading
      
      if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
        y_new <- y_init + y_change
        #print(y_new)
      } else {
        y_new <- y_init - y_change
        #print(y_new)
      }
      
      if (heading > 0 && heading <= 180) {
        x_new <- x_init + x_change
        #print(x_new)
      } else {
        x_new <- x_init - x_change
        #print(x_new)
      }
      
      y_init <- y_new
      x_init <- x_new
      
      dist <- sqrt(x_new^2 + y_new^2)
      ang.temp <- atan2(y_new, x_new)
      ang.temp <- ang.temp/0.0174532925
      if (ang.temp > 0 && ang.temp <= 90) {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      } else if (ang.temp > 90 && ang.temp <= 180) {
        ang.temp <- (360 - ang.temp)
        #print(angle)
      } else if (ang.temp < 0 && ang.temp >= -90) {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      } else {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      }
      
      df.nomad[i,1] <- x_init #lat
      df.nomad[i,2] <- y_init #long
      df.nomad[i,3] <- dist
      df.nomad[i,4] <- ang.temp
    } 
    else if (state.df[i,2] == 'Search') {
      #x_init <- nomad.df[i-1,1]
      #y_init <- nomad.df[i-1,2]
      
      random <- runif(1,0,1)
      if (random < 0.85) {
        step <- rnorm(1,450,100)
      } else {
        step <- runif(1,0,300)
      }
      
      # Uniform distribution for direction
      heading <- rnorm(1,df.nomad[i,4],60)
      
      # Determine the angle to perform trig
      angle <- 0 
      
      if (heading > 0 && heading <= 90) {
        angle <- (90 - heading)
        #print(angle)
      } else if (heading > 90 && heading <= 180) {
        angle <- (heading - 90)
        #print(angle)
      } else if (heading > 180 && heading <= 270) {
        angle <- (270 - heading)
        #print(angle)
      } else {
        angle <- (heading - 270)
        #print(angle)
      }
      
      # Use sine to determine the movement in y (latitude)
      rad_y <- angle*0.0174532925
      y_change <- sin(rad_y)*step
      y_change
      
      # Use cosine to determine the movement in x (longitude)
      rad_x <- angle*0.0174532925
      x_change <- cos(rad_x)*step
      x_change
      
      #Determine whether to add or subtract new value based on heading
      
      if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
        y_new <- y_init + y_change
        #print(y_new)
      } else {
        y_new <- y_init - y_change
        #print(y_new)
      }
      
      if (heading > 0 && heading <= 180) {
        x_new <- x_init + x_change
        #print(x_new)
      } else {
        x_new <- x_init - x_change
        #print(x_new)
      }
      
      y_init <- y_new
      x_init <- x_new
      
      dist <- sqrt(x_new^2 + y_new^2)
      ang.temp <- atan2(y_new, x_new)
      ang.temp <- ang.temp/0.0174532925
      if (ang.temp > 0 && ang.temp <= 90) {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      } else if (ang.temp > 90 && ang.temp <= 180) {
        ang.temp <- (360 - ang.temp)
        #print(angle)
      } else if (ang.temp < 0 && ang.temp >= -90) {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      } else {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      }
      
      df.nomad[i,1] <- x_init #lat
      df.nomad[i,2] <- y_init #long
      df.nomad[i,3] <- dist
      df.nomad[i,4] <- ang.temp
    }
  } else {
    
    step <- rgamma(1, shape = 2, rate = 2) * 100
    
    # Uniform distribution for direction
    heading <- runif(1,0.0001,360)
    
    # Determine the angle to perform trig
    angle <- 0 
    
    if (heading > 0 && heading <= 90) {
      angle <- (90 - heading)
      #print(angle)
    } else if (heading > 90 && heading <= 180) {
      angle <- (heading - 90)
      #print(angle)
    } else if (heading > 180 && heading <= 270) {
      angle <- (270 - heading)
      #print(angle)
    } else {
      angle <- (heading - 270)
      #print(angle)
    }
    
    # Use sine to determine the movement in y (latitude)
    rad_y <- angle*0.0174532925
    y_change <- sin(rad_y)*step
    y_change
    
    # Use cosine to determine the movement in x (longitude)
    rad_x <- angle*0.0174532925
    x_change <- cos(rad_x)*step
    x_change
    
    #Determine whether to add or subtract new value based on heading
    
    if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
      y_new <- y_init + y_change
      #print(y_new)
    } else {
      y_new <- y_init - y_change
      #print(y_new)
    }
    
    if (heading > 0 && heading <= 180) {
      x_new <- x_init + x_change
      #print(x_new)
    } else {
      x_new <- x_init - x_change
      #print(x_new)
    }
    
    y_init <- y_new
    x_init <- x_new
    
    dist <- sqrt(x_new^2 + y_new^2)
    ang.temp <- atan2(y_new, x_new)
    ang.temp <- ang.temp/0.0174532925
    if (ang.temp > 0 && ang.temp <= 90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else if (ang.temp > 90 && ang.temp <= 180) {
      ang.temp <- (360 - ang.temp)
      #print(angle)
    } else if (ang.temp < 0 && ang.temp >= -90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    }
    
    df.nomad[i,1] <- x_init #lat
    df.nomad[i,2] <- y_init #long
    df.nomad[i,3] <- dist
    df.nomad[i,4] <- ang.temp
  }
}
start.date <- '2015-01-01' 
start.time <- '00:00:00' 
interval <- 60 

increment.mins <- interval * 60 
x <- paste(start.date, start.time) 

for(i in 1:n.steps) { 
  df.nomad[i,5] <- as.POSIXct(strptime(x, "%Y-%m-%d %H:%M:%S")  + i*increment.mins)
} 

class(df.nomad[,5]) <- 'POSIXct'

plot(df.nomad[,1], df.nomad[,2], type='b')
hist(df.nomad[,3])

colnames(df.nomad) <- c("Longitude","Latitude","Dist_from_Center","Angle_from Center","Datetime")


########################################################################
############################# Territorialists ##########################
# Set starting coordinates, can be offset if spatially explicit
x_init <- 0
y_init <- 0
n.steps <- 3600
HR <- 10000
core <- 8500
dist <- 0
ang.temp <- 0

# Latitude: 111,111 meters = 1 degree
# Longitude: 111,111 meters * cos(latitude) = 1 degree

df.territorial <- data.frame(matrix(n.steps, 4))

for (i in 1:n.steps){
  
  # Function normally if distance from origin is less than home range size 
  if (dist < HR && dist < core) {
    
    step <- rgamma(1, shape = 2, rate = 2) * 750
    
    # Uniform distribution for direction
    heading <- runif(1,0.0001,360)
    
    # Determine the angle to perform trig
    angle <- 0 
    
    if (heading > 0 && heading <= 90) {
      angle <- (90 - heading)
      #print(angle)
    } else if (heading > 90 && heading <= 180) {
      angle <- (heading - 90)
      #print(angle)
    } else if (heading > 180 && heading <= 270) {
      angle <- (270 - heading)
      #print(angle)
    } else {
      angle <- (heading - 270)
      #print(angle)
    }
    
    # Use sine to determine the movement in y (latitude)
    rad_y <- angle*0.0174532925
    y_change <- sin(rad_y)*step
    y_change
    
    # Use cosine to determine the movement in x (longitude)
    rad_x <- angle*0.0174532925
    x_change <- cos(rad_x)*step
    x_change
    
    #Determine whether to add or subtract new value based on heading
    
    if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
      y_new <- y_init + y_change
      #print(y_new)
    } else {
      y_new <- y_init - y_change
      #print(y_new)
    }
    
    if (heading > 0 && heading <= 180) {
      x_new <- x_init + x_change
      #print(x_new)
    } else {
      x_new <- x_init - x_change
      #print(x_new)
    }
    
    y_init <- y_new
    x_init <- x_new
    
    dist <- sqrt(x_new^2 + y_new^2)
    ang.temp <- atan2(y_new, x_new)
    ang.temp <- ang.temp/0.0174532925
    if (ang.temp > 0 && ang.temp <= 90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else if (ang.temp > 90 && ang.temp <= 180) {
      ang.temp <- (360 - ang.temp)
      #print(angle)
    } else if (ang.temp < 0 && ang.temp >= -90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    }
    
    df.territorial[i,1] <- x_init #lat
    df.territorial[i,2] <- y_init #long
    df.territorial[i,3] <- dist
    df.territorial[i,4] <- ang.temp
    
  }
  
  # If distance is greater than home range, turn around from relative position  
  else if (dist > HR && dist > core) { 
    
    if (ang.temp < 180){
      ang.return <- 180 + ang.temp
    } else {
      ang.return <- ang.temp - 180
    }
    
    # As maintain heading towards origin if outside of core area
    
    step <- rgamma(1, shape = 2, rate = 2) * 500
    
    heading <- rnorm(1,ang.return,45)
    
    angle <- 0 
    
    if (heading > 0 && heading <= 90) {
      angle <- (90 - heading)
      #print(angle)
    } else if (heading > 90 && heading <= 180) {
      angle <- (heading - 90)
      #print(angle)
    } else if (heading > 180 && heading <= 270) {
      angle <- (270 - heading)
      #print(angle)
    } else {
      angle <- (heading - 270)
      #print(angle)
    }
    
    # Use sine to determine the movement in y (latitude)
    rad_y <- angle*0.0174532925
    y_change <- sin(rad_y)*step
    y_change
    
    # Use cosine to determine the movement in x (longitude)
    rad_x <- angle*0.0174532925
    x_change <- cos(rad_x)*step
    x_change
    
    #Determine whether to add or subtract new value based on heading
    
    if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
      y_new <- y_init + y_change
      #print(y_new)
    } else {
      y_new <- y_init - y_change
      #print(y_new)
    }
    
    if (heading > 0 && heading <= 180) {
      x_new <- x_init + x_change
      #print(x_new)
    } else {
      x_new <- x_init - x_change
      #print(x_new)
    }
    
    y_init <- y_new
    x_init <- x_new
    
    dist <- sqrt(x_new^2 + y_new^2)
    ang.temp <- atan2(y_new, x_new)
    ang.temp <- ang.temp/0.0174532925
    if (ang.temp > 0 && ang.temp <= 90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else if (ang.temp > 90 && ang.temp <= 180) {
      ang.temp <- (360 - ang.temp)
      #print(angle)
    } else if (ang.temp < 0 && ang.temp >= -90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    }
    
    df.territorial[i,1] <- x_init #lat
    df.territorial[i,2] <- y_init #long
    df.territorial[i,3] <- dist
    df.territorial[i,4] <- ang.temp
  }
  
  # If inside of the core area, return to normal functioning
  else if (dist < HR && dist > core) {if (ang.temp < 180){
    ang.return <- 180 + ang.temp
  } else {
    ang.return <- ang.temp - 180
  }
    
    # As maintain heading towards origin if outside of core area
    
    step <- rgamma(1, shape = 2, rate = 2) * 500
    
    heading <- rnorm(1,ang.return,45)
    
    angle <- 0 
    
    if (heading > 0 && heading <= 90) {
      angle <- (90 - heading)
      #print(angle)
    } else if (heading > 90 && heading <= 180) {
      angle <- (heading - 90)
      #print(angle)
    } else if (heading > 180 && heading <= 270) {
      angle <- (270 - heading)
      #print(angle)
    } else {
      angle <- (heading - 270)
      #print(angle)
    }
    
    # Use sine to determine the movement in y (latitude)
    rad_y <- angle*0.0174532925
    y_change <- sin(rad_y)*step
    y_change
    
    # Use cosine to determine the movement in x (longitude)
    rad_x <- angle*0.0174532925
    x_change <- cos(rad_x)*step
    x_change
    
    #Determine whether to add or subtract new value based on heading
    
    if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
      y_new <- y_init + y_change
      #print(y_new)
    } else {
      y_new <- y_init - y_change
      #print(y_new)
    }
    
    if (heading > 0 && heading <= 180) {
      x_new <- x_init + x_change
      #print(x_new)
    } else {
      x_new <- x_init - x_change
      #print(x_new)
    }
    
    y_init <- y_new
    x_init <- x_new
    
    dist <- sqrt(x_new^2 + y_new^2)
    ang.temp <- atan2(y_new, x_new)
    ang.temp <- ang.temp/0.0174532925
    if (ang.temp > 0 && ang.temp <= 90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else if (ang.temp > 90 && ang.temp <= 180) {
      ang.temp <- (360 - ang.temp)
      #print(angle)
    } else if (ang.temp < 0 && ang.temp >= -90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    }
    
    df.territorial[i,1] <- x_init #lat
    df.territorial[i,2] <- y_init #long
    df.territorial[i,3] <- dist
    df.territorial[i,4] <- ang.temp
  }
}

start.date <- '2015-01-01' 
start.time <- '00:00:00' 
interval <- 60 

increment.mins <- interval * 60 
x <- paste(start.date, start.time) 

for(i in 1:n.steps) { 
  df.territorial[i,5] <- as.POSIXct(strptime(x, "%Y-%m-%d %H:%M:%S")  + i*increment.mins)
} 

class(df.territorial$V5) <- 'POSIXct'

plot(df.territorial[,1], df.territorial[,2], type='b')
hist(df.territorial[,3])

colnames(df.territorial) <- c("Longitude","Latitude","Dist_from_Center","Angle_from Center","Datetime")

########################################################################
################################# Migrants #############################
# Set starting coordinates, can be offset if spatially explicit
x_init <- 0
y_init <- 0
n.steps <- 8760
dist <- 0
ang.temp <- 0

state.df <- data.frame(matrix(0,n.steps,1))
df.migrant <- data.frame(matrix(0,n.steps,4))

rand1 <- rnorm(1,3240,168)
rand2 <- rnorm(1,1104,120)
rand3 <- rnorm(1,3312,168)

reside1 <- rand1
migrate1 <- reside1 + rand2
reside2 <- migrate1 + rand3
migrate2 <- 8760 - reside2

for (i in 1:n.steps) {
  if (i < reside1) {
    state.df[i,1] <- 'Resident'
  } else if (i >= reside1 && i < migrate1) {
    state.df[i,1] <- 'Migrant1'
  } else if (i >= migrate1 && i < reside2) {
    state.df[i,1] <- 'Resident'
  } else if (i >= reside2) {
    state.df[i,1] <- 'Migrant2'
  }  
}


for (i in 1:n.steps) {
  if (i > 1) {
    if (state.df[i,1] == 'Resident') {
      
      step <- rgamma(1, shape = 2, rate = 2) * 100
      
      # Uniform distribution for direction
      heading <- runif(1,0.0001,360)
      
      # Determine the angle to perform trig
      angle <- 0 
      
      if (heading > 0 && heading <= 90) {
        angle <- (90 - heading)
        #print(angle)
      } else if (heading > 90 && heading <= 180) {
        angle <- (heading - 90)
        #print(angle)
      } else if (heading > 180 && heading <= 270) {
        angle <- (270 - heading)
        #print(angle)
      } else {
        angle <- (heading - 270)
        #print(angle)
      }
      
      # Use sine to determine the movement in y (latitude)
      rad_y <- angle*0.0174532925
      y_change <- sin(rad_y)*step
      y_change
      
      # Use cosine to determine the movement in x (longitude)
      rad_x <- angle*0.0174532925
      x_change <- cos(rad_x)*step
      x_change
      
      #Determine whether to add or subtract new value based on heading
      
      if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
        y_new <- y_init + y_change
        #print(y_new)
      } else {
        y_new <- y_init - y_change
        #print(y_new)
      }
      
      if (heading > 0 && heading <= 180) {
        x_new <- x_init + x_change
        #print(x_new)
      } else {
        x_new <- x_init - x_change
        #print(x_new)
      }
      
      y_init <- y_new
      x_init <- x_new
      
      dist <- sqrt(x_new^2 + y_new^2)
      ang.temp <- atan2(y_new, x_new)
      ang.temp <- ang.temp/0.0174532925
      if (ang.temp > 0 && ang.temp <= 90) {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      } else if (ang.temp > 90 && ang.temp <= 180) {
        ang.temp <- (360 - ang.temp)
        #print(angle)
      } else if (ang.temp < 0 && ang.temp >= -90) {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      } else {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      }
      
      df.migrant[i,1] <- x_init #lat
      df.migrant[i,2] <- y_init #long
      df.migrant[i,3] <- dist
      df.migrant[i,4] <- ang.temp
      df.migrant[i,5] <- heading
      df.migrant[i,6] <- state.df[i,1]
    } 
    
    else if (state.df[i,1] == 'Migrant1' && state.df[i-1,1] == 'Resident') {
      rand1 <- runif(1,0,1)
      if (rand1 < 0.85) {
        step <- rnorm(1,1000,100)
      } else {
        step <- runif(1,0,500)
      }
      
      # Uniform distribution for direction
      heading <- runif(1,0.0001,360)
      
      # Determine the angle to perform trig
      angle <- 0 
      
      if (heading > 0 && heading <= 90) {
        angle <- (90 - heading)
        #print(angle)
      } else if (heading > 90 && heading <= 180) {
        angle <- (heading - 90)
        #print(angle)
      } else if (heading > 180 && heading <= 270) {
        angle <- (270 - heading)
        #print(angle)
      } else {
        angle <- (heading - 270)
        #print(angle)
      }
      
      # Use sine to determine the movement in y (latitude)
      rad_y <- angle*0.0174532925
      y_change <- sin(rad_y)*step
      y_change
      
      # Use cosine to determine the movement in x (longitude)
      rad_x <- angle*0.0174532925
      x_change <- cos(rad_x)*step
      x_change
      
      #Determine whether to add or subtract new value based on heading
      
      if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
        y_new <- y_init + y_change
        #print(y_new)
      } else {
        y_new <- y_init - y_change
        #print(y_new)
      }
      
      if (heading > 0 && heading <= 180) {
        x_new <- x_init + x_change
        #print(x_new)
      } else {
        x_new <- x_init - x_change
        #print(x_new)
      }
      
      y_init <- y_new
      x_init <- x_new
      
      dist <- sqrt(x_new^2 + y_new^2)
      ang.temp <- atan2(y_new, x_new)
      ang.temp <- ang.temp/0.0174532925
      if (ang.temp > 0 && ang.temp <= 90) {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      } else if (ang.temp > 90 && ang.temp <= 180) {
        ang.temp <- (360 - ang.temp)
        #print(angle)
      } else if (ang.temp < 0 && ang.temp >= -90) {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      } else {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      }
      
      df.migrant[i,1] <- x_init #lat
      df.migrant[i,2] <- y_init #long
      df.migrant[i,3] <- dist
      df.migrant[i,4] <- ang.temp
      df.migrant[i,5] <- heading
      df.migrant[i,6] <- state.df[i,1]
    }
    
    else if (state.df[i,1] == 'Migrant1' && state.df[i-1,1] == 'Migrant1') {
      
      rand1 <- runif(1,0,1)
      if (rand1 < 0.85) {
        step <- rnorm(1,1000,100)
      } else {
        step <- runif(1,0,500)
      }
      
      # Uniform distribution for direction
      heading <- rnorm(1,df.migrant[i-1,5],1)
      
      # Determine the angle to perform trig
      angle <- 0 
      
      if (heading > 0 && heading <= 90) {
        angle <- (90 - heading)
        #print(angle)
      } else if (heading > 90 && heading <= 180) {
        angle <- (heading - 90)
        #print(angle)
      } else if (heading > 180 && heading <= 270) {
        angle <- (270 - heading)
        #print(angle)
      } else {
        angle <- (heading - 270)
        #print(angle)
      }
      
      # Use sine to determine the movement in y (latitude)
      rad_y <- angle*0.0174532925
      y_change <- sin(rad_y)*step
      y_change
      
      # Use cosine to determine the movement in x (longitude)
      rad_x <- angle*0.0174532925
      x_change <- cos(rad_x)*step
      x_change
      
      #Determine whether to add or subtract new value based on heading
      
      if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
        y_new <- y_init + y_change
        #print(y_new)
      } else {
        y_new <- y_init - y_change
        #print(y_new)
      }
      
      if (heading > 0 && heading <= 180) {
        x_new <- x_init + x_change
        #print(x_new)
      } else {
        x_new <- x_init - x_change
        #print(x_new)
      }
      
      y_init <- y_new
      x_init <- x_new
      
      dist <- sqrt(x_new^2 + y_new^2)
      ang.temp <- atan2(y_new, x_new)
      ang.temp <- ang.temp/0.0174532925
      if (ang.temp > 0 && ang.temp <= 90) {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      } else if (ang.temp > 90 && ang.temp <= 180) {
        ang.temp <- (360 - ang.temp)
        #print(angle)
      } else if (ang.temp < 0 && ang.temp >= -90) {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      } else {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      }
      
      df.migrant[i,1] <- x_init #lat
      df.migrant[i,2] <- y_init #long
      df.migrant[i,3] <- dist
      df.migrant[i,4] <- ang.temp
      df.migrant[i,5] <- heading
      df.migrant[i,6] <- state.df[i,1]
    }
    
    else if (state.df[i,1] == 'Migrant2' && state.df[i-1,1] == 'Resident') {
      rand1 <- runif(1,0,1)
      if (rand1 < 0.85) {
        step <- rnorm(1,1000,100)
      } else {
        step <- runif(1,0,500)
      }
      
      # Needs to return along the same general path
      if (ang.temp < 180){
        ang.return <- 180 + ang.temp
      } else {
        ang.return <- ang.temp - 180
      }
      
      heading <- ang.return
      
      # Determine the angle to perform trig
      angle <- 0 
      
      if (heading > 0 && heading <= 90) {
        angle <- (90 - heading)
        #print(angle)
      } else if (heading > 90 && heading <= 180) {
        angle <- (heading - 90)
        #print(angle)
      } else if (heading > 180 && heading <= 270) {
        angle <- (270 - heading)
        #print(angle)
      } else {
        angle <- (heading - 270)
        #print(angle)
      }
      
      # Use sine to determine the movement in y (latitude)
      rad_y <- angle*0.0174532925
      y_change <- sin(rad_y)*step
      y_change
      
      # Use cosine to determine the movement in x (longitude)
      rad_x <- angle*0.0174532925
      x_change <- cos(rad_x)*step
      x_change
      
      #Determine whether to add or subtract new value based on heading
      
      if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
        y_new <- y_init + y_change
        #print(y_new)
      } else {
        y_new <- y_init - y_change
        #print(y_new)
      }
      
      if (heading > 0 && heading <= 180) {
        x_new <- x_init + x_change
        #print(x_new)
      } else {
        x_new <- x_init - x_change
        #print(x_new)
      }
      
      y_init <- y_new
      x_init <- x_new
      
      dist <- sqrt(x_new^2 + y_new^2)
      ang.temp <- atan2(y_new, x_new)
      ang.temp <- ang.temp/0.0174532925
      if (ang.temp > 0 && ang.temp <= 90) {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      } else if (ang.temp > 90 && ang.temp <= 180) {
        ang.temp <- (360 - ang.temp)
        #print(angle)
      } else if (ang.temp < 0 && ang.temp >= -90) {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      } else {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      }
      
      df.migrant[i,1] <- x_init #lat
      df.migrant[i,2] <- y_init #long
      df.migrant[i,3] <- dist
      df.migrant[i,4] <- ang.temp
      df.migrant[i,5] <- heading
      df.migrant[i,6] <- state.df[i,1]
    }
    
    else if (state.df[i,1] == 'Migrant2' && state.df[i-1,1] == 'Migrant2') {
      rand1 <- runif(1,0,1)
      if (rand1 < 0.85) {
        step <- rnorm(1,1000,100)
      } else {
        step <- runif(1,0,500)
      }
      
      # Need to return along same trajectory
      if (ang.temp < 180){
        ang.return <- 180 + ang.temp
      } else {
        ang.return <- ang.temp - 180
      }
      
      heading <- rnorm(1,ang.return,5)
      
      # Determine the angle to perform trig
      angle <- 0 
      
      if (heading > 0 && heading <= 90) {
        angle <- (90 - heading)
        #print(angle)
      } else if (heading > 90 && heading <= 180) {
        angle <- (heading - 90)
        #print(angle)
      } else if (heading > 180 && heading <= 270) {
        angle <- (270 - heading)
        #print(angle)
      } else {
        angle <- (heading - 270)
        #print(angle)
      }
      
      # Use sine to determine the movement in y (latitude)
      rad_y <- angle*0.0174532925
      y_change <- sin(rad_y)*step
      y_change
      
      # Use cosine to determine the movement in x (longitude)
      rad_x <- angle*0.0174532925
      x_change <- cos(rad_x)*step
      x_change
      
      #Determine whether to add or subtract new value based on heading
      
      if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
        y_new <- y_init + y_change
        #print(y_new)
      } else {
        y_new <- y_init - y_change
        #print(y_new)
      }
      
      if (heading > 0 && heading <= 180) {
        x_new <- x_init + x_change
        #print(x_new)
      } else {
        x_new <- x_init - x_change
        #print(x_new)
      }
      
      y_init <- y_new
      x_init <- x_new
      
      dist <- sqrt(x_new^2 + y_new^2)
      ang.temp <- atan2(y_new, x_new)
      ang.temp <- ang.temp/0.0174532925
      if (ang.temp > 0 && ang.temp <= 90) {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      } else if (ang.temp > 90 && ang.temp <= 180) {
        ang.temp <- (360 - ang.temp)
        #print(angle)
      } else if (ang.temp < 0 && ang.temp >= -90) {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      } else {
        ang.temp <- (90 - ang.temp)
        #print(angle)
      }
      
      df.migrant[i,1] <- x_init #lat
      df.migrant[i,2] <- y_init #long
      df.migrant[i,3] <- dist
      df.migrant[i,4] <- ang.temp
      df.migrant[i,5] <- heading
      df.migrant[i,6] <- state.df[i,1]
    }
  }   
  else {
    
    step <- rgamma(1, shape = 2, rate = 2) * 100
    
    # Uniform distribution for direction
    heading <- runif(1,0.0001,360)
    
    # Determine the angle to perform trig
    angle <- 0 
    
    if (heading > 0 && heading <= 90) {
      angle <- (90 - heading)
      #print(angle)
    } else if (heading > 90 && heading <= 180) {
      angle <- (heading - 90)
      #print(angle)
    } else if (heading > 180 && heading <= 270) {
      angle <- (270 - heading)
      #print(angle)
    } else {
      angle <- (heading - 270)
      #print(angle)
    }
    
    # Use sine to determine the movement in y (latitude)
    rad_y <- angle*0.0174532925
    y_change <- sin(rad_y)*step
    y_change
    
    # Use cosine to determine the movement in x (longitude)
    rad_x <- angle*0.0174532925
    x_change <- cos(rad_x)*step
    x_change
    
    #Determine whether to add or subtract new value based on heading
    
    if (heading > 270 && heading <= 360 || heading > 0 && heading <= 90) {
      y_new <- y_init + y_change
      #print(y_new)
    } else {
      y_new <- y_init - y_change
      #print(y_new)
    }
    
    if (heading > 0 && heading <= 180) {
      x_new <- x_init + x_change
      #print(x_new)
    } else {
      x_new <- x_init - x_change
      #print(x_new)
    }
    
    y_init <- y_new
    x_init <- x_new
    
    dist <- sqrt(x_new^2 + y_new^2)
    ang.temp <- atan2(y_new, x_new)
    ang.temp <- ang.temp/0.0174532925
    if (ang.temp > 0 && ang.temp <= 90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else if (ang.temp > 90 && ang.temp <= 180) {
      ang.temp <- (360 - ang.temp)
      #print(angle)
    } else if (ang.temp < 0 && ang.temp >= -90) {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    } else {
      ang.temp <- (90 - ang.temp)
      #print(angle)
    }
    
    df.migrant[i,1] <- x_init #lat
    df.migrant[i,2] <- y_init #long
    df.migrant[i,3] <- dist
    df.migrant[i,4] <- ang.temp
    df.migrant[i,5] <- heading
    df.migrant[i,6] <- state.df[i,1]
  }
}

start.date <- '2015-01-01' 
start.time <- '00:00:00' 
interval <- 60 

increment.mins <- interval * 60 
x <- paste(start.date, start.time) 

for(i in 1:n.steps) { 
  df.migrant[i,7] <- as.POSIXct(strptime(x, "%Y-%m-%d %H:%M:%S")  + i*increment.mins)
} 

class(df.migrant[,7]) <- 'POSIXct'

plot(df.migrant[,1], df.migrant[,2], type='b')
hist(df.migrant[,3])

colnames(df.migrant) <- c("Longitude","Latitude","Dist_from_Center","Angle_from Center","Heading", "State", "Datetime")
