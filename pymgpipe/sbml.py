"""
SBML import and export using python-libsbml.

- The SBML importer supports all versions of SBML and the fbc package.
- The SBML exporter writes SBML L3 models.
- Annotation information is stored on the cobrapy objects
- Information from the group package is read

Parsing of fbc models was implemented as efficient as possible, whereas
(discouraged) fallback solutions are not optimized for efficiency.

Notes are only supported in a minimal way relevant for constraint-based
models. I.e., structured information from notes in the form
   <p>key: value</p>
is read into the Object.notes dictionary when reading SBML files.
On writing the Object.notes dictionary is serialized to the SBML
notes information.

Annotations are read in the Object.annotation fields.

Some SBML related issues are still open, please refer to the respective issue:
- update annotation format and support qualifiers (depends on decision
    for new annotation format; https://github.com/opencobra/cobrapy/issues/684)
- write compartment annotations and notes (depends on updated first-class
    compartments; see https://github.com/opencobra/cobrapy/issues/760)
- support compression on file handles (depends on solution for
    https://github.com/opencobra/cobrapy/issues/812)
"""

from __future__ import absolute_import

import datetime
import logging
import os
import re
from collections import defaultdict, namedtuple
from copy import deepcopy
from sys import platform

import libsbml

import cobra
from cobra.core import Gene, Group, Metabolite, Model, Reaction
from cobra.core.gene import ast2str, parse_gpr
from cobra.manipulation.validate import check_metabolite_compartment_formula
from cobra.util.solver import linear_reaction_coefficients, set_objective


try:
    from cStringIO import StringIO  # Python 2
except ImportError:
    from io import StringIO


class CobraSBMLError(Exception):
    """SBML error class."""

    pass


LOGGER = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
# Defaults and constants for writing SBML
# -----------------------------------------------------------------------------
config = cobra.Configuration()  # for default bounds
LOWER_BOUND_ID = "cobra_default_lb"
UPPER_BOUND_ID = "cobra_default_ub"
ZERO_BOUND_ID = "cobra_0_bound"

BOUND_MINUS_INF = "minus_inf"
BOUND_PLUS_INF = "plus_inf"

SBO_FBA_FRAMEWORK = "SBO:0000624"
SBO_DEFAULT_FLUX_BOUND = "SBO:0000626"
SBO_FLUX_BOUND = "SBO:0000625"
SBO_EXCHANGE_REACTION = "SBO:0000627"

LONG_SHORT_DIRECTION = {"maximize": "max", "minimize": "min"}
SHORT_LONG_DIRECTION = {"min": "minimize", "max": "maximize"}

Unit = namedtuple("Unit", ["kind", "scale", "multiplier", "exponent"])
UNITS_FLUX = (
    "mmol_per_gDW_per_hr",
    [
        Unit(kind=libsbml.UNIT_KIND_MOLE, scale=-3, multiplier=1, exponent=1),
        Unit(kind=libsbml.UNIT_KIND_GRAM, scale=0, multiplier=1, exponent=-1),
        Unit(kind=libsbml.UNIT_KIND_SECOND, scale=0, multiplier=3600, exponent=-1),
    ],
)

# -----------------------------------------------------------------------------
# Functions for id replacements (import/export)
# -----------------------------------------------------------------------------
SBML_DOT = "__SBML_DOT__"

# -----------------------------------------------------------------------------
# Precompiled note pattern
# -----------------------------------------------------------------------------
pattern_notes = re.compile(
    r"<(?P<prefix>(\w+:)?)p[^>]*>(?P<content>.*?)</(?P=prefix)p>",
    re.IGNORECASE | re.DOTALL,
)

pattern_to_sbml = re.compile(r"([^0-9_a-zA-Z])")

pattern_from_sbml = re.compile(r"__(\d+)__")


def _escape_non_alphanum(nonASCII):
    """converts a non alphanumeric character to a string representation of
    its ascii number"""
    return "__" + str(ord(nonASCII.group())) + "__"


def _number_to_chr(numberStr):
    """converts an ascii number to a character"""
    return chr(int(numberStr.group(1)))


def _clip(sid, prefix):
    """Clips a prefix from the beginning of a string if it exists."""
    return sid[len(prefix) :] if sid.startswith(prefix) else sid


def _f_gene(sid, prefix="G_"):
    """Clips gene prefix from id."""
    sid = sid.replace(SBML_DOT, ".")
    sid = pattern_from_sbml.sub(_number_to_chr, sid)
    return _clip(sid, prefix)


def _f_gene_rev(sid, prefix="G_"):
    """Adds gene prefix to id."""
    sid = pattern_to_sbml.sub(_escape_non_alphanum, sid)
    return prefix + sid.replace(".", SBML_DOT)


def _f_specie(sid, prefix="M_"):
    """Clips specie/metabolite prefix from id."""
    sid = pattern_from_sbml.sub(_number_to_chr, sid)
    return _clip(sid, prefix)


def _f_specie_rev(sid, prefix="M_"):
    """Adds specie/metabolite prefix to id."""
    sid = pattern_to_sbml.sub(_escape_non_alphanum, sid)
    return prefix + sid


def _f_reaction(sid, prefix="R_"):
    """Clips reaction prefix from id."""
    sid = pattern_from_sbml.sub(_number_to_chr, sid)
    return _clip(sid, prefix)


def _f_reaction_rev(sid, prefix="R_"):
    """Adds reaction prefix to id."""
    sid = pattern_to_sbml.sub(_escape_non_alphanum, sid)
    return prefix + sid


def _f_group(sid, prefix="G_"):
    """Clips group prefix from id."""
    sid = pattern_from_sbml.sub(_number_to_chr, sid)
    return _clip(sid, prefix)


def _f_group_rev(sid, prefix="G_"):
    """Adds group prefix to id."""
    sid = pattern_to_sbml.sub(_escape_non_alphanum, sid)
    return prefix + sid


F_GENE = "F_GENE"
F_GENE_REV = "F_GENE_REV"
F_SPECIE = "F_SPECIE"
F_SPECIE_REV = "F_SPECIE_REV"
F_REACTION = "F_REACTION"
F_REACTION_REV = "F_REACTION_REV"
F_GROUP = "F_GROUP"
F_GROUP_REV = "F_GROUP_REV"

F_REPLACE = {
    F_GENE: _f_gene,
    F_GENE_REV: _f_gene_rev,
    F_SPECIE: _f_specie,
    F_SPECIE_REV: _f_specie_rev,
    F_REACTION: _f_reaction,
    F_REACTION_REV: _f_reaction_rev,
    F_GROUP: _f_group,
    F_GROUP_REV: _f_group_rev,
}

# -----------------------------------------------------------------------------
# Write SBML
# -----------------------------------------------------------------------------
def write_sbml_model(cobra_model, filename, f_replace=F_REPLACE, **kwargs):
    """Writes cobra model to filename.

    The created model is SBML level 3 version 1 (L1V3) with
    fbc package v2 (fbc-v2).

    If the given filename ends with the suffix ".gz" (for example,
    "myfile.xml.gz"), libSBML assumes the caller wants the file to be
    written compressed in gzip format. Similarly, if the given filename
    ends with ".zip" or ".bz2", libSBML assumes the caller wants the
    file to be compressed in zip or bzip2 format (respectively). Files
    whose names lack these suffixes will be written uncompressed. Special
    considerations for the zip format: If the given filename ends with
    ".zip", the file placed in the zip archive will have the suffix
    ".xml" or ".sbml".  For example, the file in the zip archive will
    be named "test.xml" if the given filename is "test.xml.zip" or
    "test.zip". Similarly, the filename in the archive will be
    "test.sbml" if the given filename is "test.sbml.zip".

    Parameters
    ----------
    cobra_model : cobra.core.Model
        Model instance which is written to SBML
    filename : string
        path to which the model is written
    f_replace: dict of replacement functions for id replacement
    """
    doc = _model_to_sbml(cobra_model, f_replace=f_replace, **kwargs)

    if isinstance(filename, str):
        # write to path
        libsbml.writeSBMLToFile(doc, filename)

    elif hasattr(filename, "write"):
        # write to file handle
        sbml_str = libsbml.writeSBMLToString(doc)
        filename.write(sbml_str)


def _model_to_sbml(cobra_model, f_replace=None, units=True):
    """Convert Cobra model to SBMLDocument.

    Parameters
    ----------
    cobra_model : cobra.core.Model
        Cobra model instance
    f_replace : dict of replacement functions
        Replacement to apply on identifiers.
    units : boolean
        Should the FLUX_UNITS be written in the SBMLDocument.

    Returns
    -------
    libsbml.SBMLDocument
    """
    if f_replace is None:
        f_replace = {}

    sbml_ns = libsbml.SBMLNamespaces(3, 1)  # SBML L3V1
    sbml_ns.addPackageNamespace("fbc", 2)  # fbc-v2

    doc = libsbml.SBMLDocument(sbml_ns)  # type: libsbml.SBMLDocument
    doc.setPackageRequired("fbc", False)
    doc.setSBOTerm(SBO_FBA_FRAMEWORK)

    model = doc.createModel()  # type: libsbml.Model
    model_fbc = model.getPlugin("fbc")  # type: libsbml.FbcModelPlugin
    model_fbc.setStrict(True)

    if cobra_model.id is not None:
        model.setId(cobra_model.id)
        model.setMetaId("meta_" + cobra_model.id)
    else:
        model.setMetaId("meta_model")
    if cobra_model.name is not None:
        model.setName(cobra_model.name)

    # for parsing annotation corresponding to the model
    _sbase_annotations(model, cobra_model.annotation)
    # for parsing notes corresponding to the model
    _sbase_notes_dict(model, cobra_model.notes)

    # Meta information (ModelHistory) related to SBMLDocument
    if hasattr(cobra_model, "_sbml"):
        meta = cobra_model._sbml
        if "annotation" in meta:
            _sbase_annotations(doc, meta["annotation"])
        if "notes" in meta:
            _sbase_notes_dict(doc, meta["notes"])

        history = libsbml.ModelHistory()  # type: libsbml.ModelHistory
        if "created" in meta and meta["created"]:
            history.setCreatedDate(meta["created"])
        else:
            time = datetime.datetime.now()
            timestr = time.strftime("%Y-%m-%dT%H:%M:%S")
            date = libsbml.Date(timestr)
            _check(history.setCreatedDate(date), "set creation date")
            _check(history.setModifiedDate(date), "set modified date")

        if "creators" in meta:
            for cobra_creator in meta[
                "creators"
            ]:  # noqa: E501 type: libsbml.ModelCreator
                creator = libsbml.ModelCreator()
                if cobra_creator.get("familyName", None):
                    creator.setFamilyName(cobra_creator["familyName"])
                if cobra_creator.get("givenName", None):
                    creator.setGivenName(cobra_creator["givenName"])
                if cobra_creator.get("organisation", None):
                    creator.setOrganisation(cobra_creator["organisation"])
                if cobra_creator.get("email", None):
                    creator.setEmail(cobra_creator["email"])

                _check(history.addCreator(creator), "adding creator to ModelHistory.")

        # TODO: Will be implemented as part of
        #  https://github.com/opencobra/cobrapy/issues/810
        # _check(model.setModelHistory(history), 'set model history')

    # Units
    if units:
        flux_udef = model.createUnitDefinition()  # type:libsbml.UnitDefinition
        flux_udef.setId(UNITS_FLUX[0])
        for u in UNITS_FLUX[1]:
            unit = flux_udef.createUnit()  # type: libsbml.Unit
            unit.setKind(u.kind)
            unit.setExponent(u.exponent)
            unit.setScale(u.scale)
            unit.setMultiplier(u.multiplier)

    min_value = config.lower_bound
    max_value = config.upper_bound

    _create_parameter(
        model, pid=LOWER_BOUND_ID, value=min_value, sbo=SBO_DEFAULT_FLUX_BOUND
    )
    _create_parameter(
        model, pid=UPPER_BOUND_ID, value=max_value, sbo=SBO_DEFAULT_FLUX_BOUND
    )
    _create_parameter(model, pid=ZERO_BOUND_ID, value=0, sbo=SBO_DEFAULT_FLUX_BOUND)
    _create_parameter(
        model, pid=BOUND_MINUS_INF, value=-float("Inf"), sbo=SBO_DEFAULT_FLUX_BOUND
    )
    _create_parameter(
        model, pid=BOUND_PLUS_INF, value=float("Inf"), sbo=SBO_DEFAULT_FLUX_BOUND
    )

    # Compartments
    # FIXME: use first class compartment model (and write notes & annotations)
    #     (https://github.com/opencobra/cobrapy/issues/811)
    for cid, name in cobra_model.compartments.items():
        compartment = model.createCompartment()  # type: libsbml.Compartment
        compartment.setId(cid)
        compartment.setName(name)
        compartment.setConstant(True)

        # FIXME: write annotations and notes
        # _sbase_notes(c, com.notes)
        # _sbase_annotations(c, com.annotation)

    # Species
    for metabolite in cobra_model.metabolites:
        specie = model.createSpecies()  # type: libsbml.Species
        mid = metabolite.id
        if f_replace and F_SPECIE_REV in f_replace:
            mid = f_replace[F_SPECIE_REV](mid)
        specie.setId(mid)
        specie.setConstant(False)
        specie.setBoundaryCondition(False)
        specie.setHasOnlySubstanceUnits(False)
        specie.setName(metabolite.name)
        specie.setCompartment(metabolite.compartment)
        s_fbc = specie.getPlugin("fbc")  # type: libsbml.FbcSpeciesPlugin
        if metabolite.charge is not None:
            s_fbc.setCharge(metabolite.charge)
        if metabolite.formula is not None:
            s_fbc.setChemicalFormula(metabolite.formula)

        _sbase_annotations(specie, metabolite.annotation)
        _sbase_notes_dict(specie, metabolite.notes)

    # Genes
    for cobra_gene in cobra_model.genes:
        gp = model_fbc.createGeneProduct()  # type: libsbml.GeneProduct
        gid = cobra_gene.id
        if f_replace and F_GENE_REV in f_replace:
            gid = f_replace[F_GENE_REV](gid)
        gp.setId(gid)
        gname = cobra_gene.name
        if gname is None or len(gname) == 0:
            gname = gid
        gp.setName(gname)
        gp.setLabel(gid)

        _sbase_annotations(gp, cobra_gene.annotation)
        _sbase_notes_dict(gp, cobra_gene.notes)

    # Objective
    objective = model_fbc.createObjective()  # type: libsbml.Objective
    objective.setId("obj")
    objective.setType(SHORT_LONG_DIRECTION[cobra_model.objective.direction])
    model_fbc.setActiveObjectiveId("obj")

    # Reactions
    reaction_coefficients = linear_reaction_coefficients(cobra_model)
    for cobra_reaction in cobra_model.reactions:
        rid = cobra_reaction.id
        if f_replace and F_REACTION_REV in f_replace:
            rid = f_replace[F_REACTION_REV](rid)
        reaction = model.createReaction()  # type: libsbml.Reaction
        reaction.setId(rid)
        reaction.setName(cobra_reaction.name)
        reaction.setFast(False)
        reaction.setReversible((cobra_reaction.lower_bound < 0))
        _sbase_annotations(reaction, cobra_reaction.annotation)
        _sbase_notes_dict(reaction, cobra_reaction.notes)

        # stoichiometry
        for metabolite, stoichiometry in cobra_reaction._metabolites.items():
            sid = metabolite.id
            if f_replace and F_SPECIE_REV in f_replace:
                sid = f_replace[F_SPECIE_REV](sid)
            if stoichiometry < 0:
                sref = (
                    reaction.createReactant()
                )  # noqa: E501 type: libsbml.SpeciesReference
                sref.setSpecies(sid)
                sref.setStoichiometry(-stoichiometry)
                sref.setConstant(True)
            else:
                sref = (
                    reaction.createProduct()
                )  # noqa: E501 type: libsbml.SpeciesReference
                sref.setSpecies(sid)
                sref.setStoichiometry(stoichiometry)
                sref.setConstant(True)

        # bounds
        r_fbc = reaction.getPlugin("fbc")  # type: libsbml.FbcReactionPlugin
        r_fbc.setLowerFluxBound(
            _create_bound(
                model,
                cobra_reaction,
                "lower_bound",
                f_replace=f_replace,
                units=units,
                flux_udef=flux_udef,
            )
        )
        r_fbc.setUpperFluxBound(
            _create_bound(
                model,
                cobra_reaction,
                "upper_bound",
                f_replace=f_replace,
                units=units,
                flux_udef=flux_udef,
            )
        )

        # GPR
        gpr = cobra_reaction.gene_reaction_rule
        if gpr is not None and len(gpr) > 0:
            # we will parse the GPR to format them correctly for libSBML
            # avoids cases where a GPR can be parsed by cobrapy but not libSBML
            tree, gids = parse_gpr(gpr)
            # replace ids in string
            if f_replace and F_GENE_REV in f_replace:
                idmap = {gid: f_replace[F_GENE_REV](gid) for gid in gids}
                gpr_new = ast2str(tree, names=idmap)
            else:
                gpr_new = ast2str(tree)

            gpa = (
                r_fbc.createGeneProductAssociation()
            )  # noqa: E501 type: libsbml.GeneProductAssociation
            # uses ids to identify GeneProducts (True),
            # does not create GeneProducts (False)
            _check(gpa.setAssociation(gpr_new, True, False), "set gpr: " + gpr_new)

        # objective coefficients
        if reaction_coefficients.get(cobra_reaction, 0) != 0:
            flux_obj = (
                objective.createFluxObjective()
            )  # noqa: E501 type: libsbml.FluxObjective
            flux_obj.setReaction(rid)
            flux_obj.setCoefficient(cobra_reaction.objective_coefficient)

    # write groups
    if len(cobra_model.groups) > 0:
        doc.enablePackage(
            "http://www.sbml.org/sbml/level3/version1/groups/version1", "groups", True
        )
        doc.setPackageRequired("groups", False)
        model_group = model.getPlugin(
            "groups"
        )  # noqa: E501 type: libsbml.GroupsModelPlugin
        for cobra_group in cobra_model.groups:
            group = model_group.createGroup()  # type: libsbml.Group
            if f_replace and F_GROUP_REV in f_replace:
                gid = f_replace[F_GROUP_REV](cobra_group.id)
            else:
                gid = cobra_group.id
            group.setId(gid)
            group.setName(cobra_group.name)
            group.setKind(cobra_group.kind)

            _sbase_notes_dict(group, cobra_group.notes)
            _sbase_annotations(group, cobra_group.annotation)

            for cobra_member in cobra_group.members:
                member = group.createMember()  # type: libsbml.Member
                mid = cobra_member.id
                m_type = str(type(cobra_member))

                # id replacements
                if "Reaction" in m_type:
                    if f_replace and F_REACTION_REV in f_replace:
                        mid = f_replace[F_REACTION_REV](mid)
                if "Metabolite" in m_type:
                    if f_replace and F_SPECIE_REV in f_replace:
                        mid = f_replace[F_SPECIE_REV](mid)
                if "Gene" in m_type:
                    if f_replace and F_GENE_REV in f_replace:
                        mid = f_replace[F_GENE_REV](mid)

                member.setIdRef(mid)
                if cobra_member.name and len(cobra_member.name) > 0:
                    member.setName(cobra_member.name)

    return doc


def _create_bound(model, reaction, bound_type, f_replace, units=None, flux_udef=None):
    """Creates bound in model for given reaction.

    Adds the parameters for the bounds to the SBML model.

    Parameters
    ----------
    model : libsbml.Model
        SBML model instance
    reaction : cobra.core.Reaction
        Cobra reaction instance from which the bounds are read.
    bound_type : {LOWER_BOUND, UPPER_BOUND}
        Type of bound
    f_replace : dict of id replacement functions
    units : flux units

    Returns
    -------
    Id of bound parameter.
    """
    value = getattr(reaction, bound_type)
    if value == config.lower_bound:
        return LOWER_BOUND_ID
    elif value == 0:
        return ZERO_BOUND_ID
    elif value == config.upper_bound:
        return UPPER_BOUND_ID
    elif value == -float("Inf"):
        return BOUND_MINUS_INF
    elif value == float("Inf"):
        return BOUND_PLUS_INF
    else:
        # new parameter
        rid = reaction.id
        if f_replace and F_REACTION_REV in f_replace:
            rid = f_replace[F_REACTION_REV](rid)
        pid = rid + "_" + bound_type
        _create_parameter(
            model,
            pid=pid,
            value=value,
            sbo=SBO_FLUX_BOUND,
            units=units,
            flux_udef=flux_udef,
        )
        return pid


def _create_parameter(
    model, pid, value, sbo=None, constant=True, units=None, flux_udef=None
):
    """Create parameter in SBML model."""
    parameter = model.createParameter()  # type: libsbml.Parameter
    parameter.setId(pid)
    parameter.setValue(value)
    parameter.setConstant(constant)
    if sbo:
        parameter.setSBOTerm(sbo)
    if units:
        parameter.setUnits(flux_udef.getId())


def _check_required(sbase, value, attribute):
    """Get required attribute from SBase.

    Parameters
    ----------
    sbase : libsbml.SBase
    value : existing value
    attribute: name of attribute

    Returns
    -------
    attribute value (or value if already set)
    """

    if (value is None) or (value == ""):
        msg = "Required attribute '%s' cannot be found or parsed in " "'%s'." % (
            attribute,
            sbase,
        )
        if hasattr(sbase, "getId") and sbase.getId():
            msg += " with id '%s'" % sbase.getId()
        elif hasattr(sbase, "getName") and sbase.getName():
            msg += " with name '%s'" % sbase.getName()
        elif hasattr(sbase, "getMetaId") and sbase.getMetaId():
            msg += " with metaId '%s'" % sbase.getName()
        raise CobraSBMLError(msg)
    if attribute == "id":
        if not libsbml.SyntaxChecker.isValidSBMLSId(value):
            LOGGER.error("'%s' is not a valid SBML 'SId'." % value)

    return value


def _check(value, message):
    """
    Checks the libsbml return value and logs error messages.

    If 'value' is None, logs an error message constructed using
      'message' and then exits with status code 1. If 'value' is an integer,
      it assumes it is a libSBML return status code. If the code value is
      LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
      prints an error message constructed using 'message' along with text from
      libSBML explaining the meaning of the code, and exits with status code 1.

    """
    if value is None:
        LOGGER.error(
            "Error: LibSBML returned a null value trying " "to <" + message + ">."
        )
    elif type(value) is int:
        if value == libsbml.LIBSBML_OPERATION_SUCCESS:
            return
        else:
            LOGGER.error("Error encountered trying to <" + message + ">.")
            LOGGER.error(
                "LibSBML error code {}: {}".format(
                    str(value), libsbml.OperationReturnValue_toString(value).strip()
                )
            )
    else:
        return


# -----------------------------------------------------------------------------
# Notes
# -----------------------------------------------------------------------------
def _parse_notes_dict(sbase):
    """Creates dictionary of COBRA notes.

    Parameters
    ----------
    sbase : libsbml.SBase

    Returns
    -------
    dict of notes
    """
    notes = sbase.getNotesString()
    if notes and len(notes) > 0:
        notes_store = dict()
        for match in pattern_notes.finditer(notes):
            try:
                # Python 2.7 does not allow keywords for split.
                # Python 3 can have (":", maxsplit=1)
                key, value = match.group("content").split(":", 1)
            except ValueError:
                LOGGER.debug("Unexpected content format '{}'.", match.group("content"))
                continue
            notes_store[key.strip()] = value.strip()
        return {k: v for k, v in notes_store.items() if len(v) > 0}
    else:
        return {}


def _sbase_notes_dict(sbase, notes):
    """Set SBase notes based on dictionary.

    Parameters
    ----------
    sbase : libsbml.SBase
        SBML object to set notes on
    notes : notes object
        notes information from cobra object
    """
    if notes and len(notes) > 0:
        tokens = (
            ['<html xmlns = "http://www.w3.org/1999/xhtml" >']
            + ["<p>{}: {}</p>".format(k, v) for (k, v) in notes.items()]
            + ["</html>"]
        )
        _check(
            sbase.setNotes("\n".join(tokens)),
            "Setting notes on sbase: {}".format(sbase),
        )


# -----------------------------------------------------------------------------
# Annotations
# -----------------------------------------------------------------------------
"""
cobra annotations will be dictionaries of the form:
    object.annotation = {
        'provider' : [(qualifier, entity), ...]
    }
A concrete example for a metabolite would look like the following
    metabolite.annotation = {
        'chebi': [(isVersionOf, "CHEBI:17234), (is, "CHEBI:4167),],
        'kegg.compound': [(is, "C00031")]
    }
The providers are hereby MIRIAM registry keys for collections
https://www.ebi.ac.uk/miriam/main/collections
The qualifiers are biomodel qualifiers
https://co.mbine.org/standards/qualifiers

In the current stage the new annotation format is not completely supported yet.
"""

URL_IDENTIFIERS_PATTERN = re.compile(r"^https?://identifiers.org/(.+?)[:/](.+)")

URL_IDENTIFIERS_PREFIX = "https://identifiers.org"
QUALIFIER_TYPES = {
    "is": libsbml.BQB_IS,
    "hasPart": libsbml.BQB_HAS_PART,
    "isPartOf": libsbml.BQB_IS_PART_OF,
    "isVersionOf": libsbml.BQB_IS_VERSION_OF,
    "hasVersion": libsbml.BQB_HAS_VERSION,
    "isHomologTo": libsbml.BQB_IS_HOMOLOG_TO,
    "isDescribedBy": libsbml.BQB_IS_DESCRIBED_BY,
    "isEncodedBy": libsbml.BQB_IS_ENCODED_BY,
    "encodes": libsbml.BQB_ENCODES,
    "occursIn": libsbml.BQB_OCCURS_IN,
    "hasProperty": libsbml.BQB_HAS_PROPERTY,
    "isPropertyOf": libsbml.BQB_IS_PROPERTY_OF,
    "hasTaxon": libsbml.BQB_HAS_TAXON,
    "unknown": libsbml.BQB_UNKNOWN,
    "bqm_is": libsbml.BQM_IS,
    "bqm_isDescribedBy": libsbml.BQM_IS_DESCRIBED_BY,
    "bqm_isDerivedFrom": libsbml.BQM_IS_DERIVED_FROM,
    "bqm_isInstanceOf": libsbml.BQM_IS_INSTANCE_OF,
    "bqm_hasInstance": libsbml.BQM_HAS_INSTANCE,
    "bqm_unknown": libsbml.BQM_UNKNOWN,
}


def _parse_annotations(sbase):
    """Parses cobra annotations from a given SBase object.

    Annotations are dictionaries with the providers as keys.

    Parameters
    ----------
    sbase : libsbml.SBase
        SBase from which the SBML annotations are read

    Returns
    -------
    dict (annotation dictionary)

    FIXME: annotation format must be updated (this is a big collection of
          fixes) - see: https://github.com/opencobra/cobrapy/issues/684)
    """
    annotation = {}

    # SBO term
    if sbase.isSetSBOTerm():
        # FIXME: correct handling of annotations
        annotation["sbo"] = sbase.getSBOTermID()

    # RDF annotation
    cvterms = sbase.getCVTerms()
    if cvterms is None:
        return annotation

    for cvterm in cvterms:  # type: libsbml.CVTerm
        for k in range(cvterm.getNumResources()):
            # FIXME: read and store the qualifier

            uri = cvterm.getResourceURI(k)
            data = _parse_annotation_info(uri)
            if data is None:
                continue
            else:
                provider, identifier = data

            if provider in annotation:
                if isinstance(annotation[provider], str):
                    annotation[provider] = [annotation[provider]]
                # FIXME: use a list
                if identifier not in annotation[provider]:
                    annotation[provider].append(identifier)
            else:
                # FIXME: always in list
                annotation[provider] = identifier

    return annotation


def _parse_annotation_info(uri):
    """Parses provider and term from given identifiers annotation uri.

    Parameters
    ----------
    uri : str
        uri (identifiers.org url)

    Returns
    -------
    (provider, identifier) if resolvable, None otherwise
    """
    match = URL_IDENTIFIERS_PATTERN.match(uri)
    if match:
        provider, identifier = match.group(1), match.group(2)
        if provider.isupper():
            identifier = "%s:%s" % (provider, identifier)
            provider = provider.lower()
    else:
        LOGGER.warning(
            "%s does not conform to "
            "'http(s)://identifiers.org/collection/id' or"
            "'http(s)://identifiers.org/COLLECTION:id",
            uri,
        )
        return None

    return provider, identifier


def _sbase_annotations(sbase, annotation):
    """Set SBase annotations based on cobra annotations.

    Parameters
    ----------
    sbase : libsbml.SBase
        SBML object to annotate
    annotation : cobra annotation structure
        cobra object with annotation information

    FIXME: annotation format must be updated
        (https://github.com/opencobra/cobrapy/issues/684)
    """

    if not annotation or len(annotation) == 0:
        return

    # standardize annotations
    annotation_data = deepcopy(annotation)

    for key, value in annotation_data.items():
        # handling of non-string annotations (e.g. integers)
        if isinstance(value, (float, int)):
            value = str(value)
        if isinstance(value, str):
            annotation_data[key] = [("is", value)]

    for key, value in annotation_data.items():
        for idx, item in enumerate(value):
            if isinstance(item, str):
                value[idx] = ("is", item)

    # set metaId
    meta_id = "meta_{}".format(sbase.getId())
    sbase.setMetaId(meta_id)

    # rdf_items = []
    for provider, data in annotation_data.items():

        # set SBOTerm
        if provider in ["SBO", "sbo"]:
            if provider == "SBO":
                LOGGER.warning(
                    "'SBO' provider is deprecated, " "use 'sbo' provider instead"
                )
            sbo_term = data[0][1]
            _check(sbase.setSBOTerm(sbo_term), "Setting SBOTerm: {}".format(sbo_term))

            # FIXME: sbo should also be written as CVTerm
            continue

        for item in data:
            qualifier_str, entity = item[0], item[1]
            qualifier = QUALIFIER_TYPES.get(qualifier_str, None)
            if qualifier is None:
                qualifier = libsbml.BQB_IS
                LOGGER.error(
                    "Qualifier type is not supported on "
                    "annotation: '{}'".format(qualifier_str)
                )

            qualifier_type = libsbml.BIOLOGICAL_QUALIFIER
            if qualifier_str.startswith("bqm_"):
                qualifier_type = libsbml.MODEL_QUALIFIER

            cv = libsbml.CVTerm()  # type: libsbml.CVTerm
            cv.setQualifierType(qualifier_type)
            if qualifier_type == libsbml.BIOLOGICAL_QUALIFIER:
                cv.setBiologicalQualifierType(qualifier)
            elif qualifier_type == libsbml.MODEL_QUALIFIER:
                cv.setModelQualifierType(qualifier)
            else:
                raise CobraSBMLError("Unsupported qualifier: " "%s" % qualifier)
            resource = "%s/%s/%s" % (URL_IDENTIFIERS_PREFIX, provider, entity)
            cv.addResource(resource)
            _check(
                sbase.addCVTerm(cv),
                "Setting cvterm: {}, resource: {}".format(cv, resource),
            )


# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------
def validate_sbml_model(
    filename,
    check_model=True,
    internal_consistency=True,
    check_units_consistency=False,
    check_modeling_practice=False,
    **kwargs
):
    """Validate SBML model and returns the model along with a list of errors.

    Parameters
    ----------
    filename : str
        The filename (or SBML string) of the SBML model to be validated.
    internal_consistency: boolean {True, False}
        Check internal consistency.
    check_units_consistency: boolean {True, False}
        Check consistency of units.
    check_modeling_practice: boolean {True, False}
        Check modeling practise.
    check_model: boolean {True, False}
        Whether to also check some basic model properties such as reaction
        boundaries and compartment formulas.

    Returns
    -------
    (model, errors)
    model : :class:`~cobra.core.Model.Model` object
        The cobra model if the file could be read successfully or None
        otherwise.
    errors : dict
        Warnings and errors grouped by their respective types.

    Raises
    ------
    CobraSBMLError
    """
    # Errors and warnings are grouped based on their type. SBML_* types are
    # from the libsbml validator. COBRA_* types are from the cobrapy SBML
    # parser.
    keys = (
        "SBML_FATAL",
        "SBML_ERROR",
        "SBML_SCHEMA_ERROR",
        "SBML_WARNING",
        "COBRA_FATAL",
        "COBRA_ERROR",
        "COBRA_WARNING",
        "COBRA_CHECK",
    )
    errors = {key: [] for key in keys}

    # [1] libsbml validation
    doc = _get_doc_from_filename(filename)  # type: libsbml.SBMLDocument

    # set checking of units & modeling practise
    doc.setConsistencyChecks(
        libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, check_units_consistency
    )
    doc.setConsistencyChecks(
        libsbml.LIBSBML_CAT_MODELING_PRACTICE, check_modeling_practice
    )

    # check internal consistency
    if internal_consistency:
        doc.checkInternalConsistency()
    doc.checkConsistency()

    for k in range(doc.getNumErrors()):
        e = doc.getError(k)  # type: libsbml.SBMLError
        msg = _error_string(e, k=k)
        sev = e.getSeverity()
        if sev == libsbml.LIBSBML_SEV_FATAL:
            errors["SBML_FATAL"].append(msg)
        elif sev == libsbml.LIBSBML_SEV_ERROR:
            errors["SBML_ERROR"].append(msg)
        elif sev == libsbml.LIBSBML_SEV_SCHEMA_ERROR:
            errors["SBML_SCHEMA_ERROR"].append(msg)
        elif sev == libsbml.LIBSBML_SEV_WARNING:
            errors["SBML_WARNING"].append(msg)

    # [2] cobrapy validation (check that SBML can be read into model)
    # all warnings generated while loading will be logged as errors
    log_stream = StringIO()
    stream_handler = logging.StreamHandler(log_stream)
    formatter = logging.Formatter("%(levelname)s:%(message)s")
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.INFO)
    LOGGER.addHandler(stream_handler)
    LOGGER.propagate = False

    try:
        # read model and allow additional parser arguments
        model = _sbml_to_model(doc, **kwargs)
    except CobraSBMLError as e:
        errors["COBRA_ERROR"].append(str(e))
        return None, errors
    except Exception as e:
        errors["COBRA_FATAL"].append(str(e))
        return None, errors

    cobra_errors = log_stream.getvalue().split("\n")
    for cobra_error in cobra_errors:
        tokens = cobra_error.split(":")
        error_type = tokens[0]
        error_msg = ":".join(tokens[1:])

        if error_type == "WARNING":
            errors["COBRA_WARNING"].append(error_msg)
        elif error_type == "ERROR":
            errors["COBRA_ERROR"].append(error_msg)

    # remove stream handler
    LOGGER.removeHandler(stream_handler)
    LOGGER.propagate = True

    # [3] additional model tests
    if check_model:
        errors["COBRA_CHECK"].extend(check_metabolite_compartment_formula(model))

    for key in ["SBML_FATAL", "SBML_ERROR", "SBML_SCHEMA_ERROR"]:
        if len(errors[key]) > 0:
            LOGGER.error("SBML errors in validation, check error log " "for details.")
            break
    for key in ["SBML_WARNING"]:
        if len(errors[key]) > 0:
            LOGGER.error("SBML warnings in validation, check error log " "for details.")
            break
    for key in ["COBRA_FATAL", "COBRA_ERROR"]:
        if len(errors[key]) > 0:
            LOGGER.error("COBRA errors in validation, check error log " "for details.")
            break
    for key in ["COBRA_WARNING", "COBRA_CHECK"]:
        if len(errors[key]) > 0:
            LOGGER.error(
                "COBRA warnings in validation, check error log " "for details."
            )
            break

    return model, errors


def _error_string(error, k=None):
    """String representation of SBMLError.

    Parameters
    ----------
    error : libsbml.SBMLError
    k : index of error

    Returns
    -------
    string representation of error
    """
    package = error.getPackage()
    if package == "":
        package = "core"

    template = "E{} ({}): {} ({}, L{}); {}; {}"
    error_str = template.format(
        k,
        error.getSeverityAsString(),
        error.getCategoryAsString(),
        package,
        error.getLine(),
        error.getShortMessage(),
        error.getMessage(),
    )
    return error_str
