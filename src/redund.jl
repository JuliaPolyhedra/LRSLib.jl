function redund(m::RepMatrix)
    # if non-negative flag is set, non-negative constraints are not input
    # explicitly, and are not checked for redundancy

    # Step 2: Find a starting cobasis from default of specified order
    #         Lin is created if necessary to hold linearity space

    if m.status == :AtNoBasis
        getfirstbasis(m)
    end

    # Pivot to a starting dictionary
    # There may have been column redundancy
    # If so the linearity space is obtained and redundant
    # columns are removed. User can access linearity space
    # from lrs_mp_matrix Lin dimensions nredundcol x d+1

    # Step 3: Test each row of the dictionary to see if it is redundant

    # note some of these may have been changed in getting initial dictionary
    m_A = unsafe_load(m.P).m_A
    d = unsafe_load(m.P).d
    lastdv = unsafe_load(m.Q).lastdv

    # rows 0..lastdv are cost, decision variables, or linearities
    # other rows need to be tested

    redset = BitSet()
    for index in (lastdv + 1):(m_A + d)
        ineq = unsafe_load_inequality(m, index) # the input inequality number corr. to this index

        status = checkindex(m, index)
        if :redundant == checkindex(m, index)
            push!(redset, ineq)
        end
    end
    redset
end

function redundi(m::RepMatrix, ineq::Int)
    if m.status == :AtNoBasis
        getfirstbasis(m)
    end

    m_A = unsafe_load(m.P).m_A
    d = unsafe_load(m.P).d
    lastdv = unsafe_load(m.Q).lastdv

    index = lastdv + findfirst(i -> ineq == unsafe_load_inequality(m, i), (lastdv + 1):(m_A + d))

    :redundant == checkindex(m, index)
end
