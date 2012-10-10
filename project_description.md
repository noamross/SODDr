Model is

 - spatially explicit, on a lattice
 - 
   
Total spore production is the spore production per tree

Per tree-spore production is a function of species, age, and time-since infection

Probability of detection is a function of species, age, and time-since infection

In each step, there are species and age-specific probabilities for

 - Advancing to next size class
 - Mortality
 - Recruitment
 - Recovery rate from disease

In each step, disease is distributed from the plot to adjacent plots via the kernel



Main issue - what about disease in the inter-plot areas

1. Most transmission - Winter rainy season
2. Most growth (point of transition from one size class to next) - spring/early summer
3. Sampling/observation of plots - summer

Force of infection = Force from within + Force from adjacent sites + external force from outside the site.

Each location has
Tree pop by species, subdivided by age class, subdivided by disease state

Disease state is S, I x time since infection, R or D

Recruitment is of two forms: seed and resprouting.  Note that seed is really recruits per adult, because it includes other seed mortality processes, etc.

E is proportion of space unoccupied for recruitment 

$1- SUM(Pop of each species*space occupied by an individual)?

TODO for October 9, 2012
 - Start with defining parameters
 - Define the state variables
 - Write up as *.Rmd
   
Set up grid (masked to deal with non-sampled grid vertices)
Dealing with edge?  Maybe variance in external force?
Dealing with connectivity?  Maybe a matrix of connectivities between grid cells?

When using the dispersal kernel, do I aggregate spores from all species or use different interspecies transmissions?
 - Is this one of my model choice mechanisms?

Other things in the model choice heirarchy
 - Tanoak age classes
 -  
   
> Within redwood forests at each site, between 6 and 15 500-m2 circular plots were randomly established at intervals .100 m @Cobb2010

> At each study site, 30 plots (500 m^2^, i.e. 25.2 m diame- ter) were established with their location based on species composition that most closely resembled stands in the central and southern redwood forest subregions (Sawyer et al. 2000). Once stands were located at a site, a random starting point c. 25 m from the edge of a red- wood stand or trail was chosen for the first plot. Addi- tional plots were 100 m from the previous plot and were placed by randomly selecting a compass bearing. of 0, 90, 180 or 270 from the initial location, so long as the bearing kept us within the redwood forest type. Plots at a number of the sites were not contiguous, as the vegetation type often changed within a few hundred metres

Data-generating process for circular plots: Ask Latimer?

parms = c(
  #First species
  sp1.classes = 4,
  # Infection rates
    sp1.inf.1 = 0.33
    sp1.inf.2 = 0.32
    sp1.inf.3 = 0.30
    sp1.inf.4 = 0.24
  #
  )

states = c(

With relatively small plots, we have to consider the possibility of additional inputs from outside.

Should I consider the plots samples from a larger grid?

How do we deal with different inter-species kernels?
 - Generate a "force from species X" value that integreates dispersal from adjacent plots
     - Does this mean the kernel is different for different species?
     - Differences in infectivity could be due to both dispersal and susceptibility
     
Outputs - how probability of outbreak changes with input value

Question: Should arrival just be an increase in force of infection, or some other poisson process?


