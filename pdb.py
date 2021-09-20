import Bio.PDB as B
import sys
from scipy.spatial.distance import pdist


class Structure:
    def __init__(self, file):
        self.models = []
        parser = B.PDBParser()
        self.structure = parser.get_structure("struct", file)
        for model in self.structure.get_models():
            self.models.append(Model(model))

        self.chains = []
        self.residues = []
        self.atoms = []
        for model in self.models:
            for chain in model.chains:
                self.chains.append(chain.chain)
                for residuum in chain.residues:
                    self.residues.append(residuum.residuum)
                    for a in residuum.atoms:
                        self.atoms.append(a.atom)

        self.modelsNum = len(self.models)
        self.chainsNum = len(self.chains)
        self.residuesNum = len(self.residues)
        self.atomsNum = len(self.atoms)
        self.width = max(pdist([a.get_coord() for a in self.atoms]))

    def search_neighbors(self, ligand, radius, wanted):
        pool = B.NeighborSearch(self.atoms)
        return pool.search(ligand.get_coord(), radius, level=wanted)

    def getSurfaceBuriedRatio(self):
        surf = ([(key, val) for (key, val) in B.HSExposure.ExposureCN(self.models[0].model) if val < 27])
        bur = ([(key, val) for (key, val) in B.HSExposure.ExposureCN(self.models[0].model) if val >= 27])
        surfDict = {}
        polarOnSurf = 0
        for (am, val) in surf:
            if am.get_resname() in surfDict.keys():
                surfDict[am.get_resname()] += 1
            else:
                surfDict[am.get_resname()] = 1
            if isPolar(am.get_resname()):
                polarOnSurf += 1
        burDict = {}
        polarBur = 0
        for (am, val) in bur:
            if am.get_resname() in burDict.keys():
                burDict[am.get_resname()] += 1
            else:
                burDict[am.get_resname()] = 1
            if isPolar(am.get_resname()):
                polarBur += 1
        print(polarBur)
        return {"buried": {"ratio": len(bur) / (len(bur) + len(surf)),
                           "buried polar": polarBur / len(bur),
                           "distribution": burDict,
                           },
                "surface": {"ratio": len(surf) / (len(bur) + len(surf)),
                            "surface polar": polarOnSurf / len(surf),
                            "distribution": surfDict,
                }}


def isPolar(aa):
    return aa in {"SER", "THR", "CYS", "ASN", "GLN", "TYR", "ASP", "GLU", "LYS", "ARG", "HIS"}


class Model:
    def __init__(self, model):
        self.chains = []
        self.model = model
        for chain in self.model.get_chains():
            self.chains.append(Chain(chain))


class Chain:
    def __init__(self, chain):
        self.residues = []
        self.chain = chain
        for residuum in self.chain.get_residues():
            self.residues.append(Residuum(residuum))


class Residuum:
    def __init__(self, residuum):
        self.atoms = []
        self.residuum = residuum
        for atom in residuum:
            self.atoms.append(Atom(atom))


class Atom:
    def __init__(self, atom):
        self.atom = atom


def main(srgv):
    strctr = Structure("data/A2a.pdb")
    surfBurDict = strctr.getSurfaceBuriedRatio()
    ratio = surfBurDict["buried"]["ratio"]
    buriedComposition = surfBurDict["buried"]["distribution"]
    surfaceComposition = surfBurDict["surface"]["distribution"]
    surfacePolar = surfBurDict["surface"]["surface polar"]
    buriedPolar = surfBurDict["buried"]["buried polar"]
    print(surfacePolar)
    print(buriedPolar)
    print(ratio)

if __name__ == '__main__':
    main(sys.argv[1:])
