A short description of how this code works:


ProposedIgChain -> IgChain -> AbChain

With Factories:

ProposedIgChain -> IgChainFactory -> IgChain -> AbChainFactory -> AbChain


A ProposedIgChain is an identified chain through the HMMs that may or may not be an Ig.
The IgChainFactory checks the ProposedIgChain to see if it is indeed some sort of Ig molecule.
It creates an IgChain which has extra information.
The AbChainFactory checks the IgChain to see if it is an antibody (not A TCR, etc.).  If it is indeed an antibody,
it creates the AbChain which has all identification of CDRs and their clusters including numbering information for
each residue in the antibody chain.

Both Chain types have domain information to deal with ScFvs properly.  This is done through the IgDomain class.