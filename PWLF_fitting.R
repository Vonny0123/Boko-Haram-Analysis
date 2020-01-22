# With the code available here one can perform piecewise linear analysis on various subsets of the 
#attached Boko Haram data and output the full results along with summaries of break points 
#and gradients between them for different numbers of break points. There is also a function 
#to produce plots from the output to help visualise the output. I've also given some 
#examples of how to use these functions to explore subsets of data.
# The goal here is to identify key points in time where Boko Haram has increased of decreased
#their activity in different regions, against different target types, in different types of event
#and any combination thereof. 

{require(plyr)
require(formattable)
require(dplyr)
require(tidyr)
require(segmented)
require(stringr)
require(ggplot2)
require(viridisLite)
require(RColorBrewer)
require(reshape2)
require(gridExtra)} #Load required packages


my_piecewise_linear <- function(data, #the preprocessed Boko Haram Data
                                event_type_filter = NULL, #which event types? vector of strings
                                sub_event_type_filter = NULL, #which sub-event types? vector of strings
                                actor_filter = NULL, #filter for keywords in actor1 and 2
                                country_filter = NULL, #filter for countries
                                max_pieces = 3){ #max number of breakpoints
  BH <- data
  
  if(!is.null(event_type_filter)){ #Filtering for desired keywords
    BH <- BH[str_detect(BH$event_type, regex(event_type_filter, ignore_case = TRUE)),]
  }
  
  if(!is.null(country_filter)){ #Filter for desired country
    BH <- BH[str_detect(BH$country, regex(country_filter, ignore_case = TRUE)),]
  }
  
  if(!is.null(sub_event_type_filter)){ #Filter for sub event type
    BH <- BH[str_detect(BH$sub_event_type, regex(sub_event_type_filter, ignore_case = TRUE)),]
  }
  
  if(!is.null(actor_filter)){ #This searches for your keywords in the name of actors
    BH <- BH[str_detect(BH$actor1, regex(actor_filter, ignore_case = TRUE))|
               str_detect(BH$actor2, regex(actor_filter, ignore_case = TRUE)),]
  }
  
  #Now going to run piecewise linear regression for the desired numbers of break points
  
  BH$event_date <- as.Date(BH$event_date, '%d/%m/%Y') #convert to datetime type
  
  BH <- BH[, c('event_date', 'fatalities')] #For pwlf fitting need only date and fatalities
  BH$event_count <- 1 #For the cumulative sum
  BH <- BH[order(BH$event_date),] #order by date
  
  BH.agg <- aggregate(BH[,c('fatalities', 'event_count')], by=list(event_date=BH$event_date), FUN=sum) #add together counts for events with same date
  
  BH.agg <- BH.agg %>% #this adds a datapoint for each day and takes the cumulative sum
    complete(event_date = seq.Date(min(event_date), max(event_date), by="day"), fill=list(fatalities=0, event_count=0))
  BH.agg$fatalities_cum <- cumsum(BH.agg$fatalities)
  BH.agg$event_count_cum <- cumsum(BH.agg$event_count)
  
  BH.agg$index <- 1:length(BH.agg$event_date) #Add index for the regression
  
  #Define linear models for fatalities and event count versus index
  BH.lm.fatalities<-lm(fatalities_cum~index,data=BH.agg)
  BH.lm.event_count<-lm(event_count_cum~index,data=BH.agg)
  
  #Define lists to hold the resulting models and summary statistics
  pwl.models.fatalities <- vector(mode='list', length=max_pieces)
  summary.fatalities <- vector(mode='list', length=max_pieces)
  gradients.fatalities <- vector(mode='list', length=max_pieces)
  pwl.models.event_count <- vector(mode='list', length=max_pieces)
  summary.event_count <- vector(mode='list', length=max_pieces)
  gradients.event_count <- vector(mode='list', length=max_pieces)
  
  for(pieces in 1:max_pieces){
    #Do for fatalities
    o.fatalities<-segmented(BH.lm.fatalities,
                            seg.Z=~index,
                            model = TRUE,
                            npsi=pieces)
    
    pwl.models.fatalities[[pieces]] <- o.fatalities
    summary.fatalities[[pieces]] <- summary(o.fatalities)
    
    model = o.fatalities
    break_point_indices <- c(1) %>% 
      append(round(model$psi[,2])) %>%
      append(length(model$fitted.values))
    model$fitted.values[break_point_indices]
    
    gradients <- vector(length=pieces+1)
    for(piece in 1:length(gradients)){
      x0 <- break_point_indices[piece]
      x1 <- break_point_indices[piece+1]
      y0 <- model$fitted.values[break_point_indices][piece]
      y1 <- model$fitted.values[break_point_indices][piece+1]
      gradients[piece] <- (y1-y0)/(x1-x0)
    }
    gradients.fatalities[[pieces]] <- gradients
    
    #Do for event count
    o.event_count<-segmented(BH.lm.event_count,
                             seg.Z=~index,
                             model = TRUE,
                             npsi=pieces)
    
    pwl.models.event_count[[pieces]] <- o.event_count
    summary.event_count[[pieces]] <- summary(o.event_count)
    
    model = o.event_count
    break_point_indices <- c(1) %>% 
      append(round(model$psi[,2])) %>%
      append(length(model$fitted.values))
    model$fitted.values[break_point_indices]
    
    gradients <- vector(length=pieces+1)
    for(piece in 1:length(gradients)){
      x0 <- break_point_indices[piece]
      x1 <- break_point_indices[piece+1]
      y0 <- model$fitted.values[break_point_indices][piece]
      y1 <- model$fitted.values[break_point_indices][piece+1]
      gradients[piece] <- (y1-y0)/(x1-x0)
    }
    gradients.event_count[[pieces]] <- gradients
  }
  
  #Get break points for fatalities
  break_points_fatalities <- vector(mode='list', length=max_pieces)
  for(pieces in 1:max_pieces){
    a <- summary.fatalities[[pieces]]$psi
    break_points_fatalities[[pieces]] <- BH.agg$event_date[BH.agg$index %in% round(a[,2])]
  }
  
  #Get break points for event_count
  break_points_event_count <- vector(mode='list', length=max_pieces)
  for(pieces in 1:max_pieces){
    a <- summary.event_count[[pieces]]$psi
    break_points_event_count[[pieces]] <- BH.agg$event_date[BH.agg$index %in% round(a[,2])]
  }
  
  get_dates_fatalities <- function(i){
    ldply(break_points_fatalities[i], as.Date.character, format="%Y-%m-%d")
  }
  
  columns <- c('One Breakpoint', 'Two Breakpoints', 'Three Breakpoints', 'Four Breakpoints', 'Five Breakpoints')
  
  table_break_points_fatalities <- ldply(1:3, get_dates_fatalities)
  names(table_break_points_fatalities) <- NULL
  table_break_points_fatalities <- t(table_break_points_fatalities)
  table_break_points_fatalities <- as.data.frame(table_break_points_fatalities)
  names(table_break_points_fatalities) <- columns[1:max_pieces]
  null_formatter <- formatter("span", style = x ~ style(color = ifelse(is.na(x), 'white', 'black')))
  table_break_points_fatalities <- formattable(data.frame(table_break_points_fatalities), list(area(col=names(table_break_points_fatalities)) ~ null_formatter), col.names=columns[1:max_pieces])
  
  get_dates_event_count <- function(i){
    ldply(break_points_event_count[i], as.Date.character, format="%Y-%m-%d")
  }
  table_break_points_event_count <- ldply(1:3, get_dates_event_count)
  names(table_break_points_event_count) <- NULL
  table_break_points_event_count <- t(table_break_points_event_count)
  table_break_points_event_count <- as.data.frame(table_break_points_event_count)
  names(table_break_points_event_count) <- columns[1:max_pieces]
  null_formatter <- formatter("span", style = x ~ style(color = ifelse(is.na(x), 'white', 'black')))
  table_break_points_event_count <- formattable(data.frame(table_break_points_event_count), list(area(col=names(table_break_points_event_count)) ~ null_formatter), col.names=columns[1:max_pieces])
  
  
  gradients.fatalities.copy <- gradients.fatalities
  gradients.event_count.copy <- gradients.event_count
  null_formatter <- formatter("span", style = x ~ style(color = ifelse(is.na(x), 'white', 'black')))
  for(piece in 1:max_pieces){
    gradients.fatalities.copy[[piece]] <- gradients.fatalities.copy[[piece]][1:(max_pieces+1)]
    gradients.event_count.copy[[piece]] <- gradients.event_count.copy[[piece]][1:(max_pieces+1)]
  }
  df <- as.data.frame(gradients.fatalities.copy)
  names(df) <- columns[1:max_pieces]
  table.gradients.fatalities <- formattable(df, list(area(col=names(df)) ~ null_formatter), col.names=columns[1:max_pieces])
  
  df <- as.data.frame(gradients.event_count.copy)
  names(df) <- columns[1:max_pieces]
  table.gradients.event_count <- formattable(df, list(area(col=names(df)) ~ null_formatter), col.names=columns[1:max_pieces])
  
  output = list(
    'aggregated.data' = BH.agg,
    'break.points.fatalities' = break_points_fatalities,
    'break.points.event_count' = break_points_event_count,
    'table.break.points.fatalities' = table_break_points_fatalities,
    'table.break.points.event_count' = table_break_points_event_count,
    'model.fatalities' = pwl.models.fatalities,
    'model.event_count' = pwl.models.event_count,
    'gradients.fatalities' = gradients.fatalities,
    'gradients.event_count' = gradients.event_count,
    'table.gradients.fatalities' = table.gradients.fatalities,
    'table.gradients.event_count' = table.gradients.event_count
  )
  
  return(output)
                                
}

#Note: my_plotter is limited to 5 curves per graph through my choice of colours I've 
#made available. More curves shouldn't really be necessary.

my_plotter <- function(pwlf, #this is the output from my_piecewise_linear()
                      fatalities_or_events, #'fatalities' or 'events', str
                      breakpoints_to_plot=NULL, #Which numbers of breakpoints to plot, vector
                      legend.position=c(.85,.25), #where to position legend, 'none' for no legend
                      breakpoint_lines=TRUE){ #Whether to plot the lines of breakpoints or not
  
  colours <- c('#332288', '#44AA99', '#117733', '#CC6677', '#AA4499') #Chosen to be distinct non-primary colours
  
  if(fatalities_or_events=='fatalities'){
    model <- pwlf$model.fatalities
    breakpoints <- pwlf$break.points.fatalities
    data <- pwlf$aggregated.data[,c('event_date', 'fatalities_cum')]
  }else{
    model <- pwlf$model.event_count
    breakpoints <- pwlf$break.points.event_count
    data <- pwlf$aggregated.data[,c('event_date', 'event_count_cum')]
  }
  
  
  for(pieces in breakpoints_to_plot){
    newcol <- as.data.frame(predict(model[[pieces]]))
    names(newcol) <- paste(toString(pieces),'breakpoints',sep=' ')
    data <- cbind(data,newcol)
  }
  names(data)[2] <- 'True Curve'
  
  data <- melt(data, id.vars='event_date')
  
  plot <- ggplot(data=data, aes(x=event_date, y=value, color=variable))+
    scale_color_manual(values=colours)+
    geom_line(size=1.2)+#, linetype='dashed', alpha=1)+
    labs(x='', 
         y='',
         title='',
         color='legend')+
    theme(legend.position=legend.position,
          legend.background=element_blank(),
          legend.title=element_blank(),
          legend.key = element_blank())
  if(breakpoint_lines){
    color_number=2
    for(number_breakpoints in breakpoints_to_plot){
      for(breakpoint in breakpoints[[number_breakpoints]]){
        plot <- plot+geom_vline(xintercept=breakpoint, color=colours[color_number], size=0.9, linetype='dashed')
      }
      color_number <- color_number+1
    }
    
  }

  return(plot)
  
}

#######################

#First set working directory to folder containing BH
boko_haram <- read.csv('BH.csv') #load the data
boko_haram$is_urban <- as.logical(boko_haram$is_urban) #convert is_urban to logical

#Here is an example of how to use the above to compare the activity of Boko Haram in
#urban and rural areas based on analysis of the cumulative sums of fatalities and event count.
#There is a faceted plot to compare side by side and rendered tables output by the above function.

#Perform piecewise linear analysis for urban and rural data separately
x_rural <- my_piecewise_linear(boko_haram[!boko_haram$is_urban,], max_pieces=3) 
x_urban <- my_piecewise_linear(boko_haram[boko_haram$is_urban,], max_pieces=3) 

{ #block of code creating nice plot to display  cumulative sums
plot_rural_event_count <- my_plotter(x_rural,'event_count', 1:3, breakpoint_lines = TRUE, legend.position = c(0.15,0.75))
plot_rural_fatalities <- my_plotter(x_rural,'fatalities', 1:3, breakpoint_lines = TRUE, legend.position = 'none')
plot_urban_event_count <- my_plotter(x_urban,'event_count', 1:3, breakpoint_lines = TRUE, legend.position = 'none')
plot_urban_fatalities <- my_plotter(x_urban,'fatalities', 1:3, breakpoint_lines = TRUE, legend.position = 'none')

grid.arrange(arrangeGrob(plot_rural_event_count, top='Event Count', left='Rural'), 
             arrangeGrob(plot_rural_fatalities, top='Fatalities'),
             arrangeGrob(plot_urban_event_count, left='Urban'), 
             arrangeGrob(plot_urban_fatalities), nrow=2)
}

#render nice tables of break points
x_rural$table.break.points.fatalities
x_urban$table.break.points.fatalities
x_rural$table.break.points.event_count
x_urban$table.break.points.event_count

#render nice tables of gradients
x_rural$table.gradients.fatalities
x_urban$table.gradients.fatalities
x_rural$table.gradients.event_count
x_urban$table.gradients.event_count

##################

#We can filter the data too. Let's do the same as above but looking just at events of type
#'violence against civilians' occuring in Nigeria

x_rural <- my_piecewise_linear(boko_haram[!boko_haram$is_urban,],
                               max_pieces=3, 
                               event_type_filter = c('violence against civilians'),
                               country_filter = c('nigeria')) #fit for urban and rural separately
x_urban <- my_piecewise_linear(boko_haram[boko_haram$is_urban,], 
                               max_pieces=3,
                               event_type_filter = c('violence against civilians'),
                               country_filter = c('nigeria')) 

{ #block of code creating nice plot to display  cumulative sums
  plot_rural_event_count <- my_plotter(x_rural,'event_count', 1:2, breakpoint_lines = TRUE, legend.position = c(0.12,0.7))
  plot_rural_fatalities <- my_plotter(x_rural,'fatalities', 1:2, breakpoint_lines = TRUE, legend.position = 'none')
  plot_urban_event_count <- my_plotter(x_urban,'event_count', 1:2, breakpoint_lines = TRUE, legend.position = 'none')
  plot_urban_fatalities <- my_plotter(x_urban,'fatalities', 1:2, breakpoint_lines = TRUE, legend.position = 'none')
  
  grid.arrange(arrangeGrob(plot_rural_event_count, top='Event Count', left='Rural'), 
               arrangeGrob(plot_rural_fatalities, top='Fatalities'),
               arrangeGrob(plot_urban_event_count, left='Urban'), 
               arrangeGrob(plot_urban_fatalities), nrow=2)
}

#render nice tables of break points
x_rural$table.break.points.fatalities
x_urban$table.break.points.fatalities
x_rural$table.break.points.event_count
x_urban$table.break.points.event_count

#render nice tables of gradients
x_rural$table.gradients.fatalities
x_urban$table.gradients.fatalities
x_rural$table.gradients.event_count
x_urban$table.gradients.event_count

##################

x <- my_piecewise_linear(boko_haram, max_pieces=3) 

{ #block of code creating nice plot to display  cumulative sums
  plot_event_count <- my_plotter(x,'event_count', 1:3, breakpoint_lines = TRUE, legend.position = 'none')
  plot_fatalities <- my_plotter(x,'fatalities', 1:3, breakpoint_lines = TRUE, legend.position = 'none')
  
  grid.arrange(arrangeGrob(plot_event_count, top='Event Count'), 
               arrangeGrob(plot_fatalities, top='Fatalities'), nrow=1)
}
