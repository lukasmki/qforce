import networkx as nx
import numpy as np
from calkeeper import Calculation

from ..molecule import Molecule
from ..qm.qm import QM
from .creator import CustomStructureCreator, CalculationStorage


class BDECreator(CustomStructureCreator):
    def __init__(self, molecule: Molecule, config, weight: int = 0, folder=None):
        super().__init__(weight, folder)

        self.mol = molecule
        self.coords = molecule.coords
        self.atomids = molecule.atomids

        self._unique_bonds = get_unique_bonds(molecule)
        self._bde = CalculationStorage(as_dict=True)
        self._bonds = {name: CalculationStorage() for name in self._unique_bonds}
        self.charge = config.qm.charge
        if self.charge != 0:
            raise ValueError("Total charge must be zero for BDE calculations")
        self.multiplicity = config.qm.multiplicity

    def _bond(self, qm: QM, software, bond_atomids) -> list[Calculation]:
        bond_hash = "_".join(
            [
                "_".join((bond_atomids + 1).astype(str)),
                software.hash(self.charge, self.multiplicity),
            ]
        )
        folder = qm.pathways.getdir("bde_frag", bond_hash, create=True)

        # create an iterator of structures
        # E(-) = E(A-B) - E(A) - E(B) - Enb(A B) - Eb(-)
        structs = []

        G = self.mol.topo.graph.copy()
        G.remove_edges_from([tuple(bond_atomids)])

        for ifrag, fragidx in enumerate(nx.connected_components(G)):
            coords, atmids = [], []
            frag = G.subgraph(fragidx)
            for inode, node_data in frag.nodes(data=True):
                coords.append(node_data["coords"])
                atmids.append(node_data["elem"])
            coords = np.array(coords)
            atmids = np.array(atmids)

            # determine multiplicity of fragment
            print(frag.nodes(data=True))

            structs.append((ifrag, (coords, atmids, 2)))

        calcs = qm.setup_bde_calculations(folder, structs)
        for calc in calcs:
            calc: Calculation
            calc.bond_hash = bond_hash

        return calcs

    def enouts(self):
        results = []
        return results

    def gradouts(self):
        return []

    def hessouts(self):
        return []

    # def setup_pre(self, qm: QM):
    #     print("BDE Setup PRE")
    #     software = qm.get_software("software")
    #     self._bde.calculations = {
    #         name: self._bond(qm, software, atomids)
    #         for name, atomids in self._unique_bonds.items()
    #     }

    # def check_pre(self):
    #     print("BDE Check PRE")
    #     return self._check(
    #         [calc for calcs in self._bde.calculations.values() for calc in calcs]
    #     )

    # def parse_pre(self, qm: QM):
    #     print("BDE Parse PRE")
    #     results = {}
    #     for name, calculations in self._bde.calculations.items():
    #         files = [calc.check() for calc in calculations]
    #         results[name] = [qm.read_energy(file) for file in files]
    #     self._bde.results = results

    def setup_main(self, qm: QM):
        print("BDE Setup MAIN")
        print(self._bonds)

        for name, calc in self._bde.calculations.items():
            print(name, calc)
            qm_out = self._bde.results[name]
            calcs = qm.setup_energy_calculations(calc.folder, qm_out, self.atomids)
            self._bonds[name].calculations = calcs

    def check_main(self):
        print("BDE Check MAIN")
        for frag in self._bonds.values():
            pass

    def parse_main(self, qm: QM):
        print("BDE Parse MAIN")
        pass


def get_unique_bonds(mol: Molecule):
    unique_bonds = {}
    for term in mol.terms["bond"]:
        bond_type = mol.topo.edge(term.atomids[0], term.atomids[1])["vers"]
        if bond_type not in unique_bonds:
            unique_bonds[bond_type] = term.atomids
    return unique_bonds
