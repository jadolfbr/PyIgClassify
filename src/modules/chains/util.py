from .IgChainFactory import IgChainFactory
from .ProposedIgChains import ProposedIgChains

from src.modules.chains.AbChainFactory import AbChainFactory


def create_ig_and_ab_chains(fasta_path, numbering_scheme = "modified_aho", pdbaa=True):
    """
    Create IgChains and AbChains from a fasta.  This method requires that the fasta id (>id) be unique due to biopython.
    We do not use it for UserPyIgClassify as I cannot guarantee that.
    """


    proposal = ProposedIgChains(fasta_path)

    ig_creator = IgChainFactory()
    ab_creator = AbChainFactory()
    ig_chain_set = ig_creator.create_ig_vset_chain_set_from_multi_fasta(proposal, pdbaa)
    ig_chains = ig_chain_set.get_chains()
    print "# of IgChains found: "+repr(len(ig_chains))

    ab_chains = []
    for ig_chain in ig_chains:
        ab_chain = ab_creator.create_ab_vset_chain(ig_chain, numbering_scheme)
        if ab_chain:
            ab_chains.append(ab_chain)

    print "# of AbChains found: "+repr(len(ab_chains))

    return ig_chains, ab_chains