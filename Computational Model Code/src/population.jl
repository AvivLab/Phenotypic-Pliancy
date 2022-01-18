function genpop(inds::Vector{Individual{T}}) where T
# If given individuals, then create a Population
# type using these input individuals
    Population(inds,
               deepcopy(inds[1])) # set the individuals and founder for this Population type creating
end


function genpop(founder::Individual)
# If given a founder individual,
    inds1 = geninds(founder)
    genpop(inds1)
end


function genpop()
# If not given founder, then generate founder and
# copy founder to make population of N individuals,
# then run these generated individuals through
# genpop again to create a Population type with
# these individuals
    inds1 = geninds()
    genpop(inds1)
end


function update(me::Population)
    update(me.individuals)
end
