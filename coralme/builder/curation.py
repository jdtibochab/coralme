import os
import re

import logging
log = logging.getLogger(__name__)

import tqdm
bar_format = '{desc:<75}: {percentage:.1f}%|{bar:10}| {n_fmt:>5}/{total_fmt:>5} [{elapsed}<{remaining}]'

import coralme
import cobra
import pandas
import json
import copy

from cobra.core.dictlist import DictList

cur_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(cur_dir, 'column_format.json'), 'r') as f:
    column_format = json.load(f)

class CurationList(DictList):
    """Stores CurationInfo instances in a cobra DictList object.

    This class stores the generated CurationInfo instances from the
    provided manual curation files.
    """
    def save(self):
        for i in self:
            i.save()

class CurationInfo(object):
    """CurationInfo class for handling manual curation files.

    This class loads manual curation files and stores them
    as properties in an Organism instance.

    As a general rule, adding information in input files will
    prevent the MEBuilder instance from updating from homology.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.
    """
    def __init__(self,
                 id,
                 org,
                 config = {},
                 file = "",
                 name = ""):

        logging.warning("Loading {}".format(name))
        self.id = id
        self.name = name
        if file:
            self.file = file
        else:
            self.file = id + ".txt"
        self.directory = org.directory
        self.org = org
        self.config = config
        self.data = self.load()
        self.org.__setattr__(id,copy.deepcopy(self.data))
        self.sep = config["sep"]

    def _modify_from_load(self):
        """Convert manual curation file into a coralME dataset"""
        return self.data
    def _modify_for_save(self):
        """Convert coralME dataset into a manual curation file"""
        # Modify this function to save the file.
        return None
    def _modify_for_create(self,df):
        """Modify dataset to create file when not provided"""
        return df

    def read(self):
        """Read manual curation file"""
        if self.config["pathtype"] == 'absolute':
            self.filepath = self.file
        else:
            self.filepath = self.directory + self.file
        if os.path.isfile(self.filepath):
            return pandas.read_csv(self.filepath,
                                   index_col=0,
                                   sep=self.config["sep"],
                                   comment='#',
                                   skip_blank_lines=True).fillna("")
        return None

    def load(self):
        """Load and convert manual curation file into a coralME dataset"""
        self.data = self.read()
        if self.data is None:
            self.org.curation_notes['org._get_manual_curation'].append({
                'msg':'No {} file found'.format(self.id),
                'importance':'low',
                'to_do':'Fill in {}'.format(self.filepath)
            })
            self.data = self._modify_for_create(self.config["create_file"])
            self.data.to_csv(self.filepath,sep=self.config["sep"])
        self.data = self._modify_from_load()
        return self.data
    def save(self):
        """Save complemented dataset from Organism for user reference"""
        mod = self._modify_for_save()
        if mod is None:
            return
        out_dir = self.directory + "reference_files/"
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        mod[self.columns[1:]].to_csv(out_dir + self.file,sep=self.sep)
    @property
    def columns(self):
        """Default columns that are coralME-compliant"""
        return column_format[self.file]
    @property
    def org_data(self):
        """Final dataset stored in Organism instance"""
        return self.org.__getattribute__(self.id)

    def _repr_html_(self) -> str:
            """Generate html representation of reaction.

            Returns
            -------
            str
                HTML representation of the reaction.
            """
            id = cobra.util.util.format_long_string(str(self.id), 500)
            file = cobra.util.util.format_long_string(str(self.file), 500)
            name = cobra.util.util.format_long_string(str(self.name), 500)
            directory = cobra.util.util.format_long_string(str(self.directory), 500)
            data = cobra.util.util.format_long_string(str(self.data), 500)
            org_data = cobra.util.util.format_long_string(str(self.org_data), 500)

            return f"""
            <table>
                <tr><td><strong>Identifier</strong></td><td>{id}</td></tr>
                <tr><td><strong>File</strong></td><td>{file}</td></tr>
                <tr><td><strong>Name</strong></td><td>{name}</td></tr>
                <tr><td><strong>Directory</strong></td><td>{directory}</td></tr>
            </table>
        """

class ReactionCorrections(CurationInfo):
    """Reads manual input to modify reactions in the M-model.

    This class creates the property "reaction_corrections" from
    the manual inputs in reaction_corrections.txt in an
    instance of Organism.

    Input here will modify reactions at the M-model stage
    before ME-model building.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    reaction_corrections.txt :
        reaction_id,name,gene_reaction_rule,reaction,notes
        COBALT2tpp,cobalt transport in via permease (no H+),BSU24740,cobalt2_e --> cobalt2_c,No notes
        ...
    """
    def __init__(self,
                 org,
                 id = "reaction_corrections",
                 config={},
                 file="reaction_corrections.txt",
                 name="Reaction corrections"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : ',',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_from_load(self):
        return self.data.T.to_dict()


class MetaboliteCorrections(CurationInfo):
    """Reads manual input to modify metabolites in the M-model.

    This class creates the property "metabolite_corrections" from
    the manual inputs in metabolite_corrections.txt in an
    instance of Organism.

    Input here will modify metabolites at the M-model stage
    before ME-model building.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    metabolite_corrections.txt :
        metabolite_id,name,formula,notes
        glc__D_e,Glucose,C6H12O6,"Added from database"
        ...
    """
    def __init__(self,
                 org,
                 id = "metabolite_corrections",
                 config={},
                 file="metabolite_corrections.txt",
                 name="Metabolite corrections"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : ',',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_from_load(self):
        return self.data.T.to_dict()

class ProteinLocation(CurationInfo):
    """Reads manual input to add protein locations.

    This class creates the property "protein_location" from
    the manual inputs in peptide_compartment_and_pathways.txt in an
    instance of Organism.

    Input here will modify protein locations, and translocation
    pathways in the ME-model.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    peptide_compartment_and_pathways.txt :
        Complex	Complex_compartment	Protein	Protein_compartment	translocase_pathway
        BSU02690-MONOMER	Plasma_Membrane	BSU02690()	Plasma_Membrane	s
        ...
    """
    def __init__(self,
                 org,
                 id = "protein_location",
                 config={},
                 file="peptide_compartment_and_pathways.txt",
                 name="Protein location"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_save(self):
        return self.org_data

class TranslocationMultipliers(CurationInfo):
    """Reads manual input to define translocation multipliers.

    This class creates the property "translocation_multipliers" from
    the manual inputs in translocation_multipliers.txt in an
    instance of Organism.

    Input here will modify how many pores are required for
    the translocation of a protein.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    translocation_multipliers.txt :
        Gene,YidC_MONOMER,TatE_MONOMER,TatA_MONOMER
        b1855,2.0,0.0,0.0
        ...
    """
    def __init__(self,
                 org,
                 id = "translocation_multipliers",
                 config={},
                 file="translocation_multipliers.txt",
                 name="Translocation multipliers"):
        if not file:
            file = id + ".txt"
        create_file = pandas.DataFrame()
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : ',',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_from_load(self):
        return self.data.to_dict()

class LipoproteinPrecursors(CurationInfo):
    """Reads manual input to add lipoprotein precursors.

    This class creates the property "lipoprotein_precursors" from
    the manual inputs in lipoprotein_precursors.txt in an
    instance of Organism.

    Input here will add lipoprotein precursors.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    lipoprotein_precursors.txt :
        name,gene
        AcrA,b0463
        ...
    """
    def __init__(self,
                 org,
                 id = "lipoprotein_precursors",
                 config={},
                 file="lipoprotein_precursors.txt",
                 name="Lipoprotein precursors"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : ',',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_from_load(self):
        return self.data.to_dict()[self.columns[1]]
    def _modify_for_save(self):
        return pandas.DataFrame.from_dict(
            {self.columns[1]:self.org_data}).rename_axis(self.columns[0])

class CleavedMethionine(CurationInfo):
    """Reads manual input to mark proteins that undergo
    N-terminal methionine cleavage.

    This class creates the property "cleaved_methionine" from
    the manual inputs in cleaved_methionine.txt in an
    instance of Organism.

    Input here will mark proteins for N-terminal methionine
    cleavage in the ME-model.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    cleaved_methionine.txt :
        cleaved_methionine_genes
        b4154
        ...
    """
    def __init__(self,
                 org,
                 id = "cleaved_methionine",
                 config={},
                 file="cleaved_methionine.txt",
                 name="Proteins with N-terminal methionine cleavage"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_from_load(self):
        return self.data.index.to_list()
    def _modify_for_save(self):
        return pandas.DataFrame.from_dict(
            {self.columns[0]:self.org_data}).set_index(self.columns[0])

class ManualComplexes(CurationInfo):
    """Reads manual input to modify or add complexes.

    This class creates the property "manual_complexes" from
    the manual inputs in protein_corrections.txt in an
    instance of Organism.

    Input here will add, modify complexes in the ME-model,
    as well as add, modify their modifications. You can
    add a complex modification ID in the replace column,
    which will remove that modified complex and replace
    it with your manually added one.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    This example adds SufBCD with their subunits and
    modifications, while removing SufBCD_mod_2fe3s(1).

    protein_corrections.txt :
        complex_id,name,genes,mod,replace
        SufBCD,SufBC2D Fe-S cluster scaffold complex,BSU32670(1) AND BSU32710(2) AND BSU32700(1),2fe2s(1),SufBCD_mod_2fe3s(1)
        ...
    """
    def __init__(self,
                 org,
                 id = "manual_complexes",
                 config={},
                 file="protein_corrections.txt",
                 name="Protein corrections"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : ',',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    @property
    def org_data(self):
        # Manual complexes must be retrieved from complexes_df
        df = self.org.complexes_df
        flag = (df["source"].str.contains(self.org.m_model.id) | df["source"].str.contains("Inferred"))
        new_complexes = df[flag]
        new_complexes = new_complexes.drop("source",axis=1)
        new_complexes["mod"] = [''] * new_complexes.shape[0]
        new_complexes["replace"] = [''] * new_complexes.shape[0]

        df = self.org.protein_mod
        new_mods = pandas.DataFrame(columns = new_complexes.columns)
        flag = (df["Source"].str.contains("Homology") & df["Source"].notna())
        protein_mod = df[flag]
        for cplx_id,row in protein_mod.iterrows():
            df = pandas.DataFrame.from_dict({
                    row['Core_enzyme']:{
                        "mod":row["Modifications"]
                    }
                })
            new_mods = pandas.concat([new_mods,df.T],axis=0)
        return pandas.concat([self.data,new_complexes,new_mods.fillna("")],axis=0).rename_axis("complex_id")

    def _modify_for_save(self):
        return self.org_data


class Sigmas(CurationInfo):
    """Reads manual input to modify or add sigma factors.

    This class creates the property "sigmas" from
    the manual inputs in sigma_factors.txt in an
    instance of Organism.

    Input here will mark proteins for N-terminal methionine
    cleavage in the ME-model.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    sigma_factors.txt :
        sigma,complex,genes,name
        RpoH_mono,RNAP_32H,b3461,"RNA polymerase, sigma 32 (sigma H) factor"
        ...
    """
    def __init__(self,
                 org,
                 id = "sigmas",
                 config={},
                 file="sigma_factors.txt",
                 name="Sigma factors"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : ',',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_save(self):
        return self.org_data

class RhoIndependent(CurationInfo):
    """Reads manual input to define genes with rho independent
    termination.

    This class creates the property "rho_independent" from
    the manual inputs in rho_independent.txt in an
    instance of Organism.

    Input here will mark genes with rho independent transcription
    termination.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    rho_independent.txt :
        id
        b0344
        ...
    """
    def __init__(self,
                 org,
                 id = "rho_independent",
                 config={},
                 file="rho_independent.txt",
                 name="Genes with rho-independent termination"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_from_load(self):
        return self.data.index.to_list()
    def _modify_for_save(self):
        return pandas.DataFrame.from_dict(
            {self.columns[0]:self.org_data}).set_index(self.columns[0])

class RNADegradosome(CurationInfo):
    """Reads manual input to add RNA degradosome composition.

    This class creates the property "rna_degradosome" from
    the manual inputs in rna_degradosome.txt in an
    instance of Organism.

    Input here will define the composition of the RNA degradosome.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    rna_degradosome.txt :
        enzymes
        Eno_dim_mod_mg2(4)
        ...
    """
    def __init__(self,
                 org,
                 id = "rna_degradosome",
                 config={},
                 file="rna_degradosome.txt",
                 name="RNA degradosome composition"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)

    def _modify_from_load(self):
        return {"rna_degradosome" : {"enzymes" : self.data.index.to_list()}}
    def _modify_for_save(self):
        l = self.org_data["rna_degradosome"]["enzymes"]
        return pandas.DataFrame.from_dict(
            {self.columns[0]:l}).set_index(self.columns[0])

class RNAModificationMachinery(CurationInfo):
    """Reads manual input to add RNA modification machinery.

    This class creates the property "rna_modification_df" from
    the manual inputs in rna_modification.txt in an
    instance of Organism.

    Input here will define enzymes that perform RNA modifications
    for either rRNA or tRNA in the ME-model.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    rna_modification.txt :
        modification	positions	type	enzymes	source
        D	16,17,20,20A,21	tRNA	DusB_mono
        ...
    """
    def __init__(self,
                 org,
                 id = "rna_modification_df",
                 config={},
                 file="rna_modification.txt",
                 name="RNA Modification machinery"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)

    def _modify_from_load(self):
        return self.data.astype(str)
    def _modify_for_save(self):
        return self.org_data

class RNAModificationTargets(CurationInfo):
    """Reads manual input to add RNA modification targets.

    This class creates the property "rna_modification_targets" from
    the manual inputs in post_transcriptional_modification_of_RNA.txt in an
    instance of Organism.

    Input here will define RNA genes that undergo modifications.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    post_transcriptional_modification_of_RNA.txt :
        bnum	position	modification
        b0202	20A	D
        ...
    """
    def __init__(self,
                 org,
                 id = "rna_modification_targets",
                 config={},
                 file="post_transcriptional_modification_of_RNA.txt",
                 name="RNA modification targets"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_save(self):
        return self.org_data

class EnzymeReactionAssociation(CurationInfo):
    """Reads manual input to specify enzyme-reaction associations.

    This class creates the property "enz_rxn_assoc_df" from
    the manual inputs in enzyme_reaction_association.txt in an
    instance of Organism.

    Input here will create the association between enzymes and
    reactions in the ME-model.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    enzyme_reaction_association.txt :
        Reaction	Complexes
        ADNt2pp	NUPG-MONOMER OR NUPC-MONOMER
        ...
    """
    def __init__(self,
                 org,
                 id = "enz_rxn_assoc_df",
                 config={},
                 file="enzyme_reaction_association.txt",
                 name="Enzyme-reaction associations"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_save(self):
        return self.org_data

class MEMetabolites(CurationInfo):
    """Reads manual input to replace metabolites in the M-model.

    This class creates the property "me_mets" from
    the manual inputs in me_metabolites.txt in an
    instance of Organism.

    Input here will mark metabolites in the M-model for replacement
    with their corrected E-matrix component.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    me_metabolites.txt :
        id	me_id	name	formula	compartment	type
        sufbcd_c	CPLX0-1341		SufBCD complex		REPLACE
        ...
    """
    def __init__(self,
                 org,
                 id = "me_mets",
                 config={},
                 file="me_metabolites.txt",
                 name="Metabolites to substitute from M-model"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_save(self):
        return self.org_data

class SubreactionMatrix(CurationInfo):
    """Reads manual input to add subreactions.

    This class creates the property "subreaction_matrix" from
    the manual inputs in subreaction_matrix.txt in an
    instance of Organism.

    Input here will define subreactions in the ME-model.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    subreaction_matrix.txt :
        Reaction	Metabolites	Stoichiometry
        mod_acetyl_c	accoa_c	-1.0
        mod_acetyl_c	coa_c	+1.0
        ...
    """
    def __init__(self,
                 org,
                 id = "subreaction_matrix",
                 config={},
                 file="subreaction_matrix.txt",
                 name="Matrix of subreaction stoichiometries"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_save(self):
        return self.org_data

class ReactionMatrix(CurationInfo):
    """Reads manual input to add reactions to the ME-model.

    This class creates the property "reaction_matrix" from
    the manual inputs in reaction_matrix.txt in an
    instance of Organism.

    Input here will define reactions directly in the
    ME-model. Definitions here will be added to the ME-model
    after processing the M-model into the ME-model.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    reaction_matrix.txt :
        Reaction	Metabolites	Stoichiometry
        Cs_cyto_import	cs_p	-1.0
        Cs_cyto_import	h_c	1.0
        Cs_cyto_import	cs_c	1.0
        Cs_cyto_import	h_p	-1.0
        ...
    """
    def __init__(self,
                 org,
                 id = "reaction_matrix",
                 config={},
                 file="reaction_matrix.txt",
                 name="Matrix of reaction stoichiometries"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)

class OrphanSpontReactions(CurationInfo):
    """Reads manual input to add reactions to the ME-model.

    This class creates the property "orphan_and_spont_reactions" from
    the manual inputs in orphan_and_spont_reactions.txt in an
    instance of Organism.

    Input here will mark reactions as orphan or spontaneous. Orphan
    reactions will be associated with CPLX_dummy, and spontaneous ones
    will not require enzymes for flux.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    orphan_and_spont_reactions.txt :
        name	description	is_reversible	is_spontaneous	subsystems
        CODH_Fe_loading	Loading of Fe	false	true
        ...
    """
    def __init__(self,
                 org,
                 id = "orphan_and_spont_reactions",
                 config={},
                 file="orphan_and_spont_reactions.txt",
                 name="Orphan and spontaneous reactions"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)

class SubsystemClassification(CurationInfo):
    """Reads manual input to classify subsystems for Keff estimation.

    This class creates the property "subsystem_classification" from
    the manual inputs in subsystem_classification.txt in an
    instance of Organism.

    Input here will classify subsystems in umbrella classifications which
    are then used to set a median Keff and correct it with the
    complex SASA.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    subsystem_classification.txt :
    subsystem	central_CE	central_AFN	intermediate	secondary	other
    S_Amino_acids_and_related_molecules	0	1	0	0	0
    ...
    """
    def __init__(self,
                 org,
                 id = "subsystem_classification",
                 config={},
                 file="subsystem_classification.txt",
                 name="Classification of subsystems"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_from_load(self):
        subsystems = set(r.subsystem for r in self.org.m_model.reactions if r.subsystem)
        df = self.data
        d = {}
        for c in df.columns:
            for s in df[df[c] == 1].index:
                d[s] = c
        return d

class TranslocationPathways(CurationInfo):
    """Reads manual input to define translocation pathways.

    This class creates the property "translocation_pathways" from
    the manual inputs in translocation_pathways.txt in an
    instance of Organism.

    Input here will define translocation pathways and their
    machinery.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    translocation_pathways.txt :
        pathway	enzyme
        sec	BSU27650-MONOMER
        sec	BSU35300-MONOMER
        sec	secYEG
        ...
    """
    def __init__(self,
                 org,
                 id = "translocation_pathways",
                 config={},
                 file="translocation_pathways.txt",
                 name="Translocation pathways"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = {}
        for p in df.index.unique():
            d[p] = {}
            pdf = df.loc[[p]]
            d[p]["enzymes"] = pdf["enzyme"].to_list()
        return d
    def _modify_for_save(self):
        df = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        for k,v in self.org_data.items():
            for i in v["enzymes"]:
                df1 = pandas.DataFrame.from_dict({k:{self.columns[1]:i}}).T.rename_axis(self.columns[0])
                df = pandas.concat([df,df1],axis=0)
        return df

class LipidModifications(CurationInfo):
    """Reads manual input to define lipid modification machinery.

    This class creates the property "lipid_modifications" from
    the manual inputs in lipid_modifications.txt in an
    instance of Organism.

    Input here will define enzymes that perform lipid
    modifications.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    lipid_modifications.txt :
        lipid_mod	enzymes
        pg_pe_160	Lgt_MONOMER,LspA_MONOMER
        ...
    """
    def __init__(self,
                 org,
                 id = "lipid_modifications",
                 config={},
                 file="lipid_modifications.txt",
                 name="Lipid modifications"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_from_load(self):
        df = self.data
        return {k:v.split(',') for k,v in df['enzymes'].to_dict().items()}
    def _modify_for_save(self):
        d = {k:','.join(v) for k,v in self.org_data.items()}
        return pandas.DataFrame.from_dict({self.columns[1]:d}).rename_axis(self.columns[0])

class StableRNAs(CurationInfo):
    """ Defines stable RNAs
    """
    # Not necessary anymore, but kept here for reference.
    def __init__(self,
                 org,
                 id = "stable_RNAs",
                 config={},
                 file="",
                 name="Stable RNAs"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_from_load(self):
        return self.data.index.to_list()


class RibosomeStoich(CurationInfo):
    """Reads manual input to define ribosome composition.

    This class creates the property "ribosome_stoich" from
    the manual inputs in ribosomal_proteins.txt in an
    instance of Organism.

    Input here will define the composition of the ribosome.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    ribosomal_proteins.txt :
        subunits	proteins
        30S	RpsD_mono,...
        50S	generic_23s_rRNAs,generic_5s_rRNAs,RplA_mono,...
    """
    def __init__(self,
                 org,
                 id = "ribosome_stoich",
                 config={},
                 file="ribosomal_proteins.txt",
                 name="Ribosomal proteins"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            {"proteins": {"30S": "generic_16s_rRNAs",
                          "50S": "generic_5s_rRNAs,generic_23s_rRNAs"}}).rename_axis('subunit')
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_from_load(self):
        df = self.data
        ribosome_stoich = copy.deepcopy(coralme.builder.dictionaries.ribosome_stoich)
        for s, row in df.iterrows():
            proteins = row["proteins"]
            if proteins:
                proteins = proteins.split(",")
                for p in proteins:
                    if s == "30S":
                        ribosome_stoich["30_S_assembly"]["stoich"][p] = 1
                    elif s == "50S":
                        ribosome_stoich["50_S_assembly"]["stoich"][p] = 1
        return ribosome_stoich
    def _modify_for_save(self):
        d = {k.split("_S_assembly")[0]+"S":','.join(list(v["stoich"].keys())) for k,v in self.org_data.items() if "assembly" in k}
        return pandas.DataFrame.from_dict({self.columns[1]:d}).rename_axis(self.columns[0])

class RibosomeSubreactions(CurationInfo):
    """Reads manual input to define ribosome subreactions.

    This class creates the property "ribosome_subreactions" from
    the manual inputs in ribosome_subreactions.txt in an
    instance of Organism.

    Input here will define enzymes that perform a ribosome
    subreaction.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    ribosome_subreactions.txt :
        subreaction	enzyme
        gtp_bound_30S_assembly_factor_phase1	BSU16650-MONOMER
        ...
    """
    def __init__(self,
                 org,
                 id = "ribosome_subreactions",
                 config={},
                 file="ribosome_subreactions.txt",
                 name="Ribosomal subreactions"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(coralme.builder.dictionaries.ribosome_subreactions.copy()).T.rename_axis('subreaction')\
                .drop("stoich",axis=1).drop("num_mods",axis=1)
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_from_load(self):
        return self.data.T.to_dict()
    def _modify_for_save(self):
        return pandas.DataFrame.from_dict(self.org_data).T.rename_axis(self.columns[0])

class GenericDict(CurationInfo):
    """Reads manual input to define generics.

    This class creates the property "generic_dict" from
    the manual inputs in generic_dict.txt in an
    instance of Organism.

    Input here will define generics.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    generic_dict.txt :
        generic_component	enzymes
        generic_16Sm4Cm1402	RsmH_mono,RsmI_mono
        ...
    """
    def __init__(self,
                 org,
                 id = "generic_dict",
                 config={},
                 file="generic_dict.txt",
                 name="Dictionary of generic complexes"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(coralme.builder.dictionaries.generics.copy()).T.rename_axis('generic_component')
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            d[k] = {}
            if v['enzymes']:
                d[k]['enzymes'] = v['enzymes'].split(",")
            else:
                d[k]['enzymes'] = []
        return d
    def _modify_for_save(self):
        d = {k:','.join(v["enzymes"]) for k,v in self.org_data.items()}
        return pandas.DataFrame.from_dict({self.columns[1]:d}).rename_axis(self.columns[0])

class AminoacidtRNASynthetase(CurationInfo):
    """Reads manual input to define amino acid tRNA ligases.

    This class creates the property "amino_acid_trna_synthetase" from
    the manual inputs in amino_acid_trna_synthetase.txt in an
    instance of Organism.

    Input here will define amino acid tRNA ligases.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    amino_acid_trna_synthetase.txt :
        amino_acid	enzyme
        ala__L_c	Ala_RS_tetra_mod_zn2(4)
        ...
    """
    def __init__(self,
                 org,
                 id = "amino_acid_trna_synthetase",
                 config={},
                 file="amino_acid_trna_synthetase.txt",
                 name="Amino acid to tRNA synthetase associations"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            {"enzyme": coralme.builder.dictionaries.amino_acid_trna_synthetase.copy()}).rename_axis('amino_acid')
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_from_load(self):
        return self.data.to_dict()['enzyme']
    def _modify_for_create(self,df):
        return df
    def _modify_for_save(self):
        d = {self.columns[1]:self.org_data}
        return pandas.DataFrame.from_dict(d).rename_axis(self.columns[0])

class PeptideReleaseFactors(CurationInfo):
    """Reads manual input to define peptide release factors.

    This class creates the property "peptide_release_factors" from
    the manual inputs in peptide_release_factors.txt in an
    instance of Organism.

    Input here will define peptide release factors.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    peptide_release_factors.txt :
        release_factor	enzyme
        UAA	generic_RF
        ...
    """
    def __init__(self,
                 org,
                 id = "peptide_release_factors",
                 config={},
                 file="peptide_release_factors.txt",
                 name="Peptide release factors"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.translation_stop_dict.copy()).T.rename_axis('release_factor')
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzyme"] = row["enzyme"]
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            d[k] = {}
            if v['enzyme']:
                d[k]['enzyme'] = v['enzyme']
            else:
                d[k]['enzyme'] = ''
        return d
    def _modify_for_save(self):
        return pandas.DataFrame.from_dict(self.org_data).T.rename_axis(self.columns[0])

class InitiationSubreactions(CurationInfo):
    """Reads manual input to define translation initiation subreactions.

    This class creates the property "initiation_subreactions" from
    the manual inputs in initiation_subreactions.txt in an
    instance of Organism.

    Input here will define translation initiation subreactions and their
    machinery.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    initiation_subreactions.txt :
        subreaction	enzymes
        Translation_initiation_factor_InfA	InfA_mono
        ...
    """
    def __init__(self,
                 org,
                 id = "initiation_subreactions",
                 config={},
                 file="initiation_subreactions.txt",
                 name="Translation initiation subreactions"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.initiation_subreactions.copy()).T.rename_axis('subreaction')\
            .drop("stoich",axis=1).drop("element_contribution",axis=1)
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
        return d
    def _modify_for_save(self):
        d = {k:','.join(v["enzymes"]) for k,v in self.org_data.items()}
        return pandas.DataFrame.from_dict({self.columns[1]:d}).rename_axis(self.columns[0])

class ElongationSubreactions(CurationInfo):
    """Reads manual input to define translation elongation subreactions.

    This class creates the property "elongation_subreactions" from
    the manual inputs in elongation_subreactions.txt in an
    instance of Organism.

    Input here will define translation elongation subreactions and their
    machinery.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    elongation_subreactions.txt :
        subreaction	enzymes
        FusA_mono_elongation	FusA_mono
        ...
    """
    def __init__(self,
                 org,
                 id = "elongation_subreactions",
                 config={},
                 file="elongation_subreactions.txt",
                 name="Translation elongation subreactions"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.elongation_subreactions.copy()).T.rename_axis('subreaction')\
            .drop("stoich",axis=1)
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
        return d
    def _modify_for_save(self):
        d = {k:','.join(v["enzymes"]) for k,v in self.org_data.items()}
        return pandas.DataFrame.from_dict({self.columns[1]:d}).rename_axis(self.columns[0])

class TerminationSubreactions(CurationInfo):
    """Reads manual input to define translation termination subreactions.

    This class creates the property "termination_subreactions" from
    the manual inputs in termination_subreactions.txt in an
    instance of Organism.

    Input here will define translation termination subreactions and their
    machinery.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    termination_subreactions.txt :
        subreaction	enzymes
        PrfA_mono_mediated_termination	PrfA_mono
        ...
    """
    def __init__(self,
                 org,
                 id = "termination_subreactions",
                 config={},
                 file="termination_subreactions.txt",
                 name="Translation termination subreactions"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.termination_subreactions.copy()).T.rename_axis('subreaction')\
            .drop("stoich",axis=1).drop("element_contribution",axis=1)
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
        return d
    def _modify_for_save(self):
        d = {k:','.join(v["enzymes"]) for k,v in self.org_data.items()}
        return pandas.DataFrame.from_dict({self.columns[1]:d}).rename_axis(self.columns[0])

class SpecialtRNASubreactions(CurationInfo):
    """Reads manual input to define special tRNA subreactions.

    This class creates the property "special_trna_subreactions" from
    the manual inputs in special_trna_subreactions.txt in an
    instance of Organism.

    Input here will define special tRNA subreactions, such as
    tRNA-Sec (selenocysteine) synthesis from tRNA-Ser.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    special_trna_subreactions.txt :
        subreaction	enzymes
        PrfA_mono_mediated_termination	PrfA_mono
        ...
    """
    def __init__(self,
                 org,
                 id = "special_trna_subreactions",
                 config={},
                 file="special_trna_subreactions.txt",
                 name="Special tRNA subreactions"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.special_trna_subreactions.copy()).T.rename_axis('subreaction')\
            .drop("stoich",axis=1).drop("element_contribution",axis=1)
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
        return d
    def _modify_for_save(self):
        d = {k:','.join(v["enzymes"]) for k,v in self.org_data.items()}
        return pandas.DataFrame.from_dict({self.columns[1]:d}).rename_axis(self.columns[0])

class TranscriptionSubreactions(CurationInfo):
    """Reads manual input to define transcription subreactions.

    This class creates the property "transcription_subreactions" from
    the manual inputs in transcription_subreactions.txt in an
    instance of Organism.

    Input here will define machinery for transcription subreactions. These
    subreactions are a set of pre-defined subreactions that are used
    in ME-models.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    transcription_subreactions.txt :
        mechanism	enzymes
        Transcription_normal_rho_independent	Mfd_mono_mod_mg2(1),NusA_mono,NusG_mono,GreA_mono,GreB_mono,RpoZ_mono_mod_mg2(1)
        ...
    """
    def __init__(self,
                 org,
                 id = "transcription_subreactions",
                 config={},
                 file="transcription_subreactions.txt",
                 name="Transcription subreactions"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.transcription_subreactions.copy()).T.rename_axis('mechanism')\
            .drop("stoich",axis=1)
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
        return d
    def _modify_for_save(self):
        d = {k:','.join(v["enzymes"]) for k,v in self.org_data.items()}
        return pandas.DataFrame.from_dict({self.columns[1]:d}).rename_axis(self.columns[0])

class SpecialModifications(CurationInfo):
    """Reads manual input to define machinery for special modifications.

    This class creates the property "special_modifications" from
    the manual inputs in special_modifications.txt in an
    instance of Organism.

    Input here will define machinery for special modifications. These
    modifications are a set of pre-defined modifications that are used
    in ME-models.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    special_modifications.txt :
        modification	enzymes	stoich
        fes_transfer	CPLX0-7617,IscA_tetra,CPLX0-7824
        ...
    """
    def __init__(self,
                 org,
                 id = "special_modifications",
                 config={},
                 file="special_modifications.txt",
                 name="Special protein modifications"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(coralme.builder.dictionaries.special_modifications.copy()).T.rename_axis('modification')\
            .drop("stoich",axis=1)
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
        return d
    def _modify_for_save(self):
        d = {k:','.join(v["enzymes"]) for k,v in self.org_data.items()}
        return pandas.DataFrame.from_dict({self.columns[1]:d}).rename_axis(self.columns[0])

class ExcisionMachinery(CurationInfo):
    """Reads manual input to define machinery for excision.

    This class creates the property "excision_machinery" from
    the manual inputs in excision_machinery.txt in an
    instance of Organism.

    Input here will define machinery for excision.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    excision_machinery.txt :
        mechanism	enzymes
        rRNA_containing	RNase_E_tetra_mod_zn2(2), ...
        ...
    """
    def __init__(self,
                 org,
                 id = "excision_machinery",
                 config={},
                 file="excision_machinery.txt",
                 name="Excision machinery"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.excision_machinery.copy()).T.rename_axis('type')
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            d[k] = {}
            if v['enzymes']:
                d[k]['enzymes'] = v['enzymes'].split(",")
            else:
                d[k]['enzymes'] = []
        return d
    def _modify_for_save(self):
        d = {k:','.join(v["enzymes"]) for k,v in self.org_data.items()}
        return pandas.DataFrame.from_dict({self.columns[1]:d}).rename_axis(self.columns[0])

class FoldingDict(CurationInfo):
    """Reads manual input to define folding pathways for proteins.

    This class creates the property "folding_dict" from
    the manual inputs in folding_dict.txt in an
    instance of Organism.

    Input here will define folding pathways for proteins.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.

    Examples
    --------
    folding_dict.txt :
        mechanism	enzymes
        GroEL_dependent_folding	b0014, ...
        ...
    """
    def __init__(self,
                 org,
                 id = "folding_dict",
                 config={},
                 file="folding_dict.txt",
                 name="Protein to folding machinery associations"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.folding_dict.copy()).T.rename_axis('mechanism')
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file,
                        name = name)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            d[k] = {}
            if v['enzymes']:
                d[k]['enzymes'] = v['enzymes'].split(",")
            else:
                d[k]['enzymes'] = []
        return d
    def _modify_for_save(self):
        d = {k:','.join(v["enzymes"]) for k,v in self.org_data.items()}
        return pandas.DataFrame.from_dict({self.columns[1]:d}).rename_axis(self.columns[0])

class MEManualCuration(object):
    """MEManualCuration class for loading manual curation from files

    This class acts as a loader for all manual curation inputs. It is
    used by Organism to retrieve manual curation in a coralME-ready
    format.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.
    """

    def __init__(self,
                 org):
        self.org = org
        self.directory = self.org.directory
        self.curation_notes = self.org.curation_notes
        self.is_reference = self.org.is_reference
        self.configuration = self.org.config

    def load_manual_curation(self):
        self.org.manual_curation = coralme.builder.curation.CurationList()
        self.org.manual_curation.append(ReactionCorrections(self.org))
        self.org.manual_curation.append(MetaboliteCorrections(self.org))
        self.org.manual_curation.append(ProteinLocation(self.org))
        self.org.manual_curation.append(TranslocationMultipliers(self.org))
        self.org.manual_curation.append(LipoproteinPrecursors(self.org))
        self.org.manual_curation.append(CleavedMethionine(self.org))
        self.org.manual_curation.append(ManualComplexes(self.org))
        self.org.manual_curation.append(Sigmas(self.org))
        self.org.manual_curation.append(RhoIndependent(self.org))
        self.org.manual_curation.append(RNADegradosome(self.org))
        self.org.manual_curation.append(RNAModificationMachinery(self.org))
        self.org.manual_curation.append(RNAModificationTargets(self.org))
        self.org.manual_curation.append(EnzymeReactionAssociation(self.org))
        self.org.manual_curation.append(MEMetabolites(self.org))
        self.org.manual_curation.append(SubreactionMatrix(self.org))
        self.org.manual_curation.append(ReactionMatrix(self.org))
        self.org.manual_curation.append(OrphanSpontReactions(self.org))
        self.org.manual_curation.append(SubsystemClassification(self.org))
        self.org.manual_curation.append(TranslocationPathways(self.org))
        self.org.manual_curation.append(LipidModifications(self.org))
        #self.org.manual_curation.append(StableRNAs(self.org))
        self.org.manual_curation.append(RibosomeStoich(self.org))
        self.org.manual_curation.append(RibosomeSubreactions(self.org))
        self.org.manual_curation.append(GenericDict(self.org))
        self.org.manual_curation.append(AminoacidtRNASynthetase(self.org))
        self.org.manual_curation.append(PeptideReleaseFactors(self.org))
        self.org.manual_curation.append(InitiationSubreactions(self.org))
        self.org.manual_curation.append(ElongationSubreactions(self.org))
        self.org.manual_curation.append(TerminationSubreactions(self.org))
        self.org.manual_curation.append(TranscriptionSubreactions(self.org))
        self.org.manual_curation.append(SpecialtRNASubreactions(self.org))
        self.org.manual_curation.append(SpecialModifications(self.org))
        self.org.manual_curation.append(ExcisionMachinery(self.org))
        self.org.manual_curation.append(FoldingDict(self.org))


class MECurator(object):
    """MECurator class for integrating additional curation.

    This class integrates additional manual curation that was
    not integrated in Organism.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.
    """

    def __init__(self,
                org):
        self.org = org

    def find_issue_with_query(self,query):
        for k,v in self.org.curation_notes.items():
            coralme.builder.helper_functions.find_issue(query,v)


def _str_to_dict(d):
    regex = ":(?=[-]?\d+(?:$|\.))"
    return (
        {re.split(regex, i)[0]: float(re.split(regex, i)[1]) for i in d.split(",")}
        if d
        else {}
    )

def _dict_to_str(d):
    return ",".join(["{}:{}".format(k, v) for k, v in d.items()])












