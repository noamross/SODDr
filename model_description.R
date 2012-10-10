#' Generate the lattice

#' Populate the lattice with initial values 
#'  - Each location should have a population value for each species, age class, and infection state (S, I+time since infection)
#'  
#' Iterative process
#' 1. Transmission
#'  - Generate force of infection in each location, FROM each species
#'    - Each species has a dispersal kernel. Force FROM a particular location is the sum of individuals of that species times the spores/individual for that species, with weightings by size class.
#'    - Force TO a location is summed up across all the kernels, but still subdivided by species
#'    - Species in each location are infected at a rate equal to the force of infection from each species times their susceptibility to each species, summed across species
#'      - Stochastic component:  There is a *probability* of infection, not a fixed proportion transitioning
#'    - Trees infected previously advance in disease age by year
#'    
#'  2. Growth
#'    - Trees advance a size class with fixed probability (deterministic)
#'    - Trees die with different probabilities based on species, infection condition
#'    - Trees recover at a rate depending on their species
#'    
#'  3. Observation
#'    - Trees are observed from a portion of the location
#'    - Probability of observing tree as infected increases with time since infection
#'    
#'  4. Recruitment
#'    - New trees from seed based on number of adults
#'    - New trees from sprout depending on number of dead in (2)
#'    
#'    
#' Questions for Alan:
#'  - Sampling from each location or assume that locations are all there is?
#'    
#' Questions for Dave/Rich
#'  - Are 1-year old trees censused?  Also, when do trees recruit from seed?  I assume spring
#'  - When you count "dead", do you get all the dead?  Or just at a certain size class. Also, Can you tell from a dead tree its cause (other than previously seeing infections?).
#'  
#'  Several sites, with 6-15 locations, each.  Need to fit some parameters heirarchically by site (and possibly per plot?)