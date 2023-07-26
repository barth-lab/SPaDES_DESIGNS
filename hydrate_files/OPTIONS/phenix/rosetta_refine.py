# LIBTBX_SET_DISPATCHER_NAME phenix.rosetta_refine

from __future__ import division
from __future__ import print_function
from phenix.automation.protocols.refine import refine_result
from phenix.automation.refinement import refinement_callback
from phenix.refinement import command_line
from phenix.utilities import citations
import phenix.rosetta.refine
from mmtbx.command_line import generate_master_phil_with_inputs
from mmtbx.validation import molprobity
import mmtbx.command_line
from libtbx.str_utils import make_header, make_sub_header
from libtbx import runtime_utils
import libtbx.phil
from libtbx.utils import Sorry, create_run_directory
import time
import os
import sys

def master_phil () :
  return generate_master_phil_with_inputs(
    #enable_automatic_twin_detection=False,
    enable_pdb_interpretation_params=False,
    phil_string="""\
rosetta {
  protocol = *full hires
    .type = choice
    .short_caption = Refinement protocol
    .caption = Full High-resolution
    .help = Refinement protocol to use, corresponding to XML scripts in the \
      Rosetta distribution.  The 'hires' protocol is less aggressive and \
      faster, but has a smaller radius of convergence.
  script = None
    .type = path
    .short_caption = Rosetta XML script path
    .help = XML file in RosettaScripts syntax describing a custom refinement \
      protocol.
  rosetta_path = None
    .type = path
    .short_caption = Rosetta distribution path
    .help = Path to Rosetta distribution.  Normally this is set by the \
      PHENIX_ROSETTA_PATH environment variable.
    .style = directory
  spades_args = None
    .type = path
    .short_caption = SPaDES HYDRATE protocol flags
    .help = filename (including path) to extra flags you want to include \
      during refinement. Uses the SPaDES hydrate movers (water placement)
  database = None
    .type = path
    .short_caption = Rosetta database path
    .help = Path to the Rosetta database.  Normally this is set by the \
      ROSETTA3_DB environment variable.
    .style = directory
  include scope phenix.rosetta.refine.rosetta_refine_params_str
  ligand_params_file = None
    .type = path
    .multiple = True
}
phenix {
  post_refine = False
    .type = bool
    .short_caption = Run additional refinement in phenix.refine
    .help = Perform additional refinement (coordinates and B-factors) in \
      phenix.refine after Rosetta is done.
  strip_hydrogens = False
    .type = bool
    .short_caption = Remove hydrogen atoms
    .help = Remove hydrogen atoms after Rosetta refinement.
  cycles = 3
    .type = int
    .short_caption = Number of phenix.refine cycles
    .input_size = 50
    .style = spinner
  tls = False
    .type = bool
    .short_caption = Use TLS refinement
  phil_file = None
    .type = path
    .short_caption = phenix.refine settings
    .style = input_file file_type:phil
    .style = noauto
  waters = False
    .type = bool
    .short_caption = Update water molecules
    .help = Add water molecules in phenix.refine after Rosetta refinement.
}
output {
  prefix = None
    .type = str
    .short_caption = Output prefix
  directory_number = None
    .type = int
  verbose = False
    .type = bool
  debug = True
    .type = bool
  gui_output_dir = None
    .type = path
    .help = Output directory (for PHENIX GUI only)
    .short_caption = Output directory
    .style = output_dir
  include scope libtbx.phil.interface.tracking_params
}
runtime {
  include scope libtbx.easy_mp.parallel_phil_str
}
""")

scripts = {
  "full" : "low_resolution_refine.xml",
  "hires" : "high_resolution_refine.xml",
}

def run (args, out=None, create_dir=False) :
  if (out is None) : out = sys.stdout
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil(),
    out=out,
    process_pdb_file=True,
    create_fmodel=True,
    prefer_anomalous=False,
    create_log_buffer=True,
    usage_string="""\
phenix.rosetta_refine model.pdb data.mtz [script=refine.xml] [options]

Run Rosetta refinement with Phenix ML X-ray target.  Requires separate
Python-enabled Rosetta installation.
""")
  fmodel = cmdline.fmodel
  pdb_hierarchy = cmdline.pdb_hierarchy
  counts = pdb_hierarchy.overall_counts()
  if (len(counts.consecutive_residue_groups_with_same_resid) > 0) :
    print("ERROR: consecutive residue groups with same resid:", file=out)
    counts.show_consecutive_residue_groups_with_same_resid(out=out,
      prefix="  ")
    raise Sorry("Please renumber the residues to be distinct.")
  params = cmdline.params
  if (params.rosetta.script is not None) :
    params.rosetta.script = os.path.abspath(params.rosetta.script)
  if (len(params.rosetta.ligand_params_file) > 0) :
    params.rosetta.ligand_params_file = [ os.path.abspath(fn)
      for fn in params.rosetta.ligand_params_file ]
  out = cmdline.log
  if (params.rosetta.rosetta_path is not None) :
    os.environ["PHENIX_ROSETTA_PATH"] = params.rosetta.rosetta_path
  scripts_cmd = phenix.rosetta.refine.find_scripts_command()
  if (scripts_cmd is None) :
    raise Sorry("""
    The RosettaScripts executable could not be located.  Please set the
    environmental variable PHENIX_ROSETTA_PATH or add the appropriate
    directory to your PATH environment variable.
  """)
  # if not phenix.rosetta.refine.find_rosetta_environment_dispatcher('phenix.rosetta_refine'):
  #   raise Sorry("""
  #   Rosetta environment variable could not be located in dispatcher.
  #   Please run the command
  #     rosetta.build_phenix_interface
  #   to add to dispatcher and compile new Rosetta scripts.
  # """)
  validate_params(params)
  make_header("ROSETTA/PHENIX X-ray refinement", out=out)
  if (params.output.prefix is None) :
    params.output.prefix = os.path.splitext(os.path.basename(
      params.input.pdb.file_name[0]))[0] + "_rosetta"
  prefix = params.output.prefix
  symm_file = prefix + ".symm"
  if (create_dir) :
    dir_name = create_run_directory(prefix="rosetta",
      default_directory_number=params.output.directory_number)
    os.chdir(dir_name)
  cmdline.start_log_file(prefix + ".log")
  print("Output directory:", file=out)
  print("  %s" % os.getcwd(), file=out)
  t_start = time.time()
  make_sub_header("Setting up input files for Rosetta", out=out)
  new_model, hetatm_chains = phenix.rosetta.refine.prepare_refinement_input(
    pdb_hierarchy=pdb_hierarchy,
    crystal_symmetry=fmodel.xray_structure,
    symm_file=symm_file,
    prefix=prefix,
    generate_new_symm_file=True,
    keep_ligands=params.rosetta.keep_ligands,
    log=out)
  script_file = params.rosetta.script
  if (script_file is None) :
    script_file = phenix.rosetta.refine.find_refinement_script(
      script_name=scripts[params.rosetta.protocol],
      rosetta_path=params.rosetta.rosetta_path)
  assert (script_file is not None)
  rosetta_result = phenix.rosetta.refine.run_refinement_trials(
    fmodel=fmodel,
    pdb_hierarchy=pdb_hierarchy,
    params=params.rosetta,
    mp_params=params.runtime,
    out=out,
    symm_file=symm_file,
    prefix=prefix,
    database=params.rosetta.database,
    script_file=script_file,
    ligand_params_files=params.rosetta.ligand_params_file,
    spades_args=params.rosetta.spades_args,
    debug=params.output.debug)
  output_file = rosetta_result.best_model
  if (output_file is None) or (not os.path.isfile(output_file)) :
    raise RuntimeError("Rosetta failed to write an output file - please "+
      "see console output for details.")
  current_model = output_file
  if (len(hetatm_chains) > 0) :
    print("Merging back heteroatoms from original input file.", file=out)
    current_model = phenix.rosetta.refine.merge_atoms(
      starting_model=current_model,
      other_chains=hetatm_chains)
  phenix_args = [
    rosetta_result.data_file,
    current_model,
    "nproc=%s" % str(params.runtime.nproc),
    "output.prefix=%s_phenix" % prefix,
    "output.serial=1",
    "output.serial_format=%d",
  ] + params.input.monomers.file_name
  if (params.phenix.post_refine) :
    make_sub_header("Running phenix.refine with optimization", out=out)
    phenix_args.extend([
      "main.number_of_macro_cycles=%d" % params.phenix.cycles,
      "optimize_adp_weight=True",
      "optimize_xyz_weight=True",
      "ncs_search.enabled=True",
      "ncs.type=torsion",
    ])
    strategy = ["individual_sites",]
    if (params.rosetta.adp_type == "group") :
      strategy.append("group_adp")
    else :
      strategy.append("individual_adp")
    if (params.phenix.tls) :
      strategy.append("tls")
      phenix_args.append("tls.find_automatically=True")
    phenix_args.append("refine.strategy=%s" % "+".join(strategy))
  else :
    make_sub_header("Running phenix.refine (with null strategy)", out=out)
    phenix_args.extend([
      "main.number_of_macro_cycles=1",
      "refine.strategy=None",
    ])
  old_stdout = sys.stdout
  phenix_args.append("--overwrite")
  cb = None
  if (not params.output.verbose) :
    phenix_args.append("--quiet")
    cb = refinement_callback(out)
  refine_obj = command_line.run(
    command_name="phenix.refine",
    args=phenix_args,
    call_back_handler=cb)
  if (old_stdout is not None) :
    sys.stdout = old_stdout
  prefix = refine_obj.params.output.prefix
  refine_pdb_file = os.path.join(os.getcwd(), "%s_1.pdb" % prefix)
  refine_mtz_file = os.path.join(os.getcwd(), "%s_1.mtz" % prefix)
  assert (os.path.isfile(refine_pdb_file) and os.path.isfile(refine_mtz_file))
  validation = molprobity.molprobity(
    model=refine_obj.model,
    fmodel=refine_obj.fmodel)
  # remove unpicklable attributes
  validation.model = None
  validation.pdb_hierarchy = None
  validation.model_statistics_geometry = None
  final_model = prefix + "_001.pdb"
  mtz_file = prefix + "_001.mtz"
  os.rename(refine_pdb_file, final_model)
  os.rename(refine_mtz_file, mtz_file)
  validation.show_summary(out=out)
  make_sub_header("Final results", out=out)
  print("", file=out)
  print("Refined model: %s" % final_model, file=out)
  print("Final maps:    %s" % mtz_file, file=out)
  print("", file=out)
  t_end = time.time()
  print("Elapsed time: %.1fs" % (t_end - t_start), file=out)
  print("", file=out)
  print("Citation:", file=out)
  citations.show_citations(["rosetta_refine"], out=out, format="cell")
  return result(
    program_name="rosetta_refine",
    job_title=params.output.job_title,
    directory=os.getcwd(),
    log_file=prefix+".log",
    pdb_files=[final_model],
    map_file=mtz_file,
    data_file=mtz_file,
    other_result=validation,
    statistics={
      "R-work (start)": rosetta_result.r_work_start,
      "R-free (start)": rosetta_result.r_free_start,
      "R-work" : validation.r_work(),
      "R-free" : validation.r_free(),
      "RMS(bonds)" : validation.rms_bonds(),
      "RMS(angles)" : validation.rms_angles(),
      "MolProbity score" : validation.molprobity_score(),
    })

class result (refine_result) :
  def get_final_stats (self) :
    return [ (kn, self.get_statistic(kn)) for kn in ["R-work", "R-free",
      "RMS(bonds)", "RMS(angles)"] ]

def validate_params (params) :
  import mmtbx.command_line
  return mmtbx.command_line.validate_input_params(params)

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    os.mkdir(self.output_dir)
    os.chdir(self.output_dir)
    return run(args=list(self.args))

if (__name__ == "__main__") :
  run(sys.argv[1:], create_dir=True)
