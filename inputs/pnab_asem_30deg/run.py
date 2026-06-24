from multiprocessing import freeze_support

from openbabel import openbabel as ob
import pnab


def main():
    run = pnab.pNAB("options.yaml")
    run.run()
    file_name = "%i_%i.pdb" % (run.results[0, 0], run.results[0, 1])

    conv = ob.OBConversion()
    mol = ob.OBMol()
    conv.ReadFile(mol, file_name)

    chains = mol.Separate()
    final_mol = ob.OBMol()

    for nchain, chain in enumerate(chains):
        delete_atom = []
        num_residues_original = chain.NumResidues()

        for residue in list(ob.OBResidueIter(chain)):
            residue.SetNum(residue.GetNum() * 2 - 1)

            new_residue1 = chain.NewResidue()
            new_residue1.SetName("GLU")
            new_residue1.SetNum(residue.GetNum() + 1)
            new_residue1.SetChain("ABCDEF"[nchain])

            new_residue2 = chain.NewResidue()
            new_residue2.SetName(residue.GetName())
            new_residue2.SetNum(residue.GetNum() + 2)
            new_residue2.SetChain("ABCDEF"[nchain])

            for atom in ob.OBResidueAtomIter(residue):
                name = residue.GetAtomID(atom)

                if name.count("'") == 0:
                    continue
                if name.count("'") == 1:
                    new_residue1.AddAtom(atom)
                    new_residue1.SetAtomID(atom, name.strip().strip("'"))
                elif name.count("'") == 2 and residue.GetIdx() + 1 != num_residues_original:
                    new_residue2.AddAtom(atom)
                    new_residue2.SetAtomID(atom, name.strip().strip("'"))
                elif name.count("'") == 2:
                    delete_atom.append(atom)

        for atom in delete_atom:
            chain.DeleteAtom(atom)

        final_mol += chain

    final_mol.SetChainsPerceived()
    conv.WriteFile(final_mol, "fixed.pdb")


if __name__ == "__main__":
    freeze_support()
    main()
