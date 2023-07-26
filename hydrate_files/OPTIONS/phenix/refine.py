
from __future__ import division
from __future__ import print_function
import libtbx.introspection
from libtbx.str_utils import make_header, format_value
from libtbx.utils import Sorry, null_out, multi_out
from libtbx import easy_run
from libtbx import adopt_init_args, Auto
import iotbx.pdb
import libtbx.load_env
import multiprocessing
import time
import os
import sys

def adjust_path_if_necessary (rosetta_path) :
  while True :
    if (os.path.basename(rosetta_path) in ["bin", "source", "main"]) :
      rosetta_path = os.path.dirname(rosetta_path)
    else :
      break
  return rosetta_path

def find_score_jd2_command(**kwds):
  return find_rosetta_command("score_jd2.", **kwds)

def find_relax_command (**kwds):
  return find_rosetta_command("relax.", **kwds)

def find_scripts_command (**kwds):
  return find_rosetta_command("rosetta_scripts.python.", **kwds)

def find_rosetta_environment_dispatcher(cmd):
  paths = os.environ.get('PATH', '')
  if paths.find(':')==-1: return False
  for path in paths.split(':'):
    if cmd in os.listdir(path):
      break
  f = open(os.path.join(path,cmd))
  lines = f.read()
  del f
  if lines.find('PHENIX_ROSETTA_PATH')==-1:
    return False
  return True

def find_rosetta_command (prefix, rosetta_path=None) :
  def search_path (path) :
    full_cmd = None
    for node in os.listdir(path) :
      if node.startswith(prefix) and (not "default" in node) :
        if (full_cmd is None) :
          full_cmd = os.path.join(path, node)
        elif (full_cmd.endswith("debug")) and (node.endswith("release")) :
          full_cmd = os.path.join(path, node)
    return full_cmd
  rosetta_cmd = None
  if (rosetta_cmd is None) :
    rosetta_bin = None
    if (rosetta_path is not None) :
      rosetta_path = adjust_path_if_necessary(rosetta_path)
      rosetta_bin = os.path.join(rosetta_path, "main", "source", "bin")
    elif ("ROSETTA_BIN" in os.environ) :
      rosetta_bin = os.environ["ROSETTA_BIN"]
    elif ("PHENIX_ROSETTA_PATH" in os.environ) :
      rosetta_path = adjust_path_if_necessary(os.environ["PHENIX_ROSETTA_PATH"])
      rosetta_bin = os.path.join(rosetta_path, "main", "source", "bin")
    if (rosetta_bin is not None) :
      rosetta_cmd = search_path(rosetta_bin)
  if (rosetta_cmd is None) :
    for path in os.environ.get("PATH", "").split(":") :
      if (not os.path.isdir(path)) :
        continue
      rosetta_cmd = search_path(path)
      if (rosetta_cmd is not None) :
        break
  return rosetta_cmd

def find_refinement_script (script_name, rosetta_path=None) :
  if (rosetta_path is None) :
    rosetta_path = os.environ.get("PHENIX_ROSETTA_PATH", None)
  if (rosetta_path is None) :
    raise Sorry("Can't find Rosetta installation - you must define this "+
      "using the environment variable PHENIX_ROSETTA_PATH, or the "+
      "parameter rosetta_path.")
  script_dir = os.path.join(rosetta_path, "main", "source", "src", "apps",
    "public", "crystal_refinement")
  if (not os.path.isdir(script_dir)) :
    raise Sorry(("Can't find directory containing refinement scripts "+
      "(expected %s).  Make "+
      "sure you have specified the Rosetta path correctly.") % script_dir)
  script_file = os.path.join(script_dir, script_name)
  if os.path.isfile(script_file) :
    return script_file
  else :
    raise Sorry("Can't find the script %s in the Rosetta installation." %
      script_name)

def find_symm_command () :
  symm_cmd = None
  relax_cmd = find_relax_command()
  assert (relax_cmd is not None)
  base_dir = os.path.dirname(os.path.dirname(relax_cmd))
  print(base_dir)
  script_file = os.path.join(base_dir, "src", "apps", "public", "symmetry",
    "make_symmdef_file.pl")
  return script_file

def run_symm_command (pdb_in, out_file):
  cmd = find_symm_command()
  assert (cmd is not None)
  args = ["/usr/bin/perl", cmd, "-m", "cryst", "-p", pdb_in]
  output = easy_run.fully_buffered(args)#.raise_if_errors()
  assert (len(output.stdout_lines) > 0)
  f = open(out_file, "w")
  for line in output.stdout_lines :
    print(line, file=f)
  f.close()

def get_process_info () :
  hostname = os.uname()[1].split(".")[0]
  pid = os.getpid()
  return "%s_%d" % (hostname, pid)


rosetta_score_params_str = """
  score_weights = *score12_full standard score13_env_hb
    .type = choice(multi=False)
"""

rosetta_relax_params_str = """
  %s
  optimize_hydrogens = True
    .type = bool
  n_cycles = 5
    .type = int
""" % rosetta_score_params_str

rosetta_refine_params_str = """
  %s
  number_of_models = 5
    .type = int(value_min=1)
    .optional = False
    .short_caption = Number of Rosetta models
  density_sampling = *fast thorough
    .type = choice
  ignore_unrecognized = False
    .type = bool
    .short_caption = Ignore unrecognized molecules
  keep_ligands = False
    .type = bool
    .short_caption = Don't remove ligands prior to Rosetta refinement
  adp_type = *individual group
    .type = choice
    .short_caption = B-factor refinement type
  map_type = *Auto 2mFo-DFc mFo dm prime_and_switch feature_enhanced composite_omit
    .type = choice
    .short_caption = Real-space refinement map type
    .help = This controls the procedure for generating a map for Rosetta \
      real-space refinement.  All map types other than 'mFo' will essentially \
      be similar to 2mFo-DFc, but with different processing steps.

""" % rosetta_score_params_str

# FIXME not currently used?
common_refine_params_str = """
#  twin_law = None
#    .type = str
#  target = *ml mlhl lsq
#    .type = choice
"""

def add_hydrogens (file_name, out=None, debug=True) :
  if (out is None) :
    out = sys.stdout
  assert os.path.isfile(file_name)
  from iotbx.file_reader import any_file
  f = any_file(file_name, force_type="pdb")
  f.assert_file_type("pdb")
  cryst = f.file_object.crystallographic_section()
  score_out = os.path.splitext(file_name)[0] + "_0001.pdb"
  args = ["-in:file:s %s" % file_name,
          "-out:output",
          "-out:pdb",
          "-out:file:o %s" % score_out,
          "-no_optH false",
          "-overwrite"]
  relax_cmd = find_relax_command()
  assert (relax_cmd is not None)
  path = os.path.dirname(relax_cmd)
  score_cmd = os.path.join(os.path.dirname(relax_cmd),
    "score." + os.path.basename(relax_cmd).split(".")[1])
  if (not os.path.exists(score_cmd)) :
    raise RuntimeError("Can't find %s in Rosetta installation!" % score_cmd)
  if (debug) :
    print("%s %s" % (score_cmd, " ".join(args)), file=out)
  score = easy_run.fully_buffered("%s %s" % (score_cmd, " ".join(args)))
  if (debug) :
    print("\n".join(score.stdout_lines), file=out)
  if (not os.path.isfile(score_out)) :
    raise RuntimeError("Adding hydrogens failed - dumping Rosetta stderr "+
      "output below.\n%s" % "\n".join(score.stderr_lines))
  pdb_lines = []
  pdb_lines.extend(cryst)
  for line in open(score_out, "r").readlines() :
    if (line.startswith("ATOM")) :
      pdb_lines.append(line.strip())
  pdb_lines.append("END")
  os.remove(score_out)
  pdb_out = os.path.splitext(file_name)[0] + "_rosetta.pdb"
  open(pdb_out, "w").write("\n".join(pdb_lines))
  return pdb_out

# XXX obsolete
def relax_simple (pdb_hierarchy, params, prefix="rosetta", out=None) :
  assert (not None in [pdb_hierarchy, params])
  if (out is None) :
    out = sys.stdout
  relax_cmd = find_relax_command()
  assert (relax_cmd is not None)
  make_header("ROSETTA fastrelax", out=out)
  tmp_file = "rosetta_in_%s.pdb" % get_process_info()
  pdb_out = open(tmp_file, "w")
  pdb_out.write(pdb_hierarchy.as_pdb_string())
  pdb_out.close()
  args = ["-in:file:s %s" % tmp_file,
          "-ex1", "-ex2aro",
          "-score:weights %s" % params.score_weights,
          "-relax:jump_move true",
          "-relax:fast",
          "-relax:default_repeats %d" % params.n_cycles,
          "-no_optH %s" % str(not params.optimize_hydrogens).lower(),
          "-out:file:o %s" % prefix,
          "-overwrite"]
  relax = easy_run.call("%s %s" % (relax_cmd, " ".join(args)))

def prepare_refinement_input (
    pdb_hierarchy,
    crystal_symmetry,
    symm_file,
    prefix,
    generate_new_symm_file=True,
    keep_ligands=False,
    log=None) :
  from iotbx import file_reader
  if (os.path.isfile(symm_file)) and (generate_new_symm_file) :
    os.remove(symm_file)
  assert (crystal_symmetry is not None)
  models = pdb_hierarchy.models()
  if (len(models) > 1) :
    raise Sorry("Multi-MODEL PDB files not supported.")
  for chain in models[0].chains() :
    n_confs = len(chain.conformers())
    if (n_confs > 1) :
      raise Sorry("Alternate conformations not supported by Rosetta.")
  rename_waters(pdb_hierarchy)
  hetatm_chains = []
  if (not keep_ligands) :
    hetatm_chains = remove_heteroatoms(pdb_hierarchy, log)
  new_model = prefix + "_in.pdb"
  f = open(new_model, "w")
  f.write(pdb_hierarchy.as_pdb_string(crystal_symmetry=crystal_symmetry))
  f.close()
  if (generate_new_symm_file) :
    run_symm_command(
      pdb_in=new_model,
      out_file=symm_file)
  if (not os.path.isfile(symm_file)) :
    raise RuntimeError("%s does not exist!" % symm_file)
  return new_model, hetatm_chains

def merge_atoms (
    starting_model,
    other_chains,
    crystal_symmetry=None) :
  from iotbx import file_reader
  f = file_reader.any_file(starting_model, force_type="pdb")
  if (crystal_symmetry is None) :
    crystal_symmetry = f.file_object.crystal_symmetry()
  hierarchy = f.file_object.hierarchy
  model = hierarchy.models()[0]
  for chain in other_chains :
    model.append_chain(chain)
  move_waters(hierarchy)
  new_model = os.path.splitext(starting_model)[0] + "_combined.pdb"
  f = open(new_model, "w")
  cryst1 = iotbx.pdb.format_cryst1_record(crystal_symmetry)
  f.write(cryst1 + "\n")
  f.write(hierarchy.as_pdb_string())
  f.close()
  return new_model

class run_refinement_trials (object) :
  def __init__ (self,
                fmodel,
                pdb_hierarchy,
                params,
                mp_params,
                out,
                symm_file,
                prefix,
                database=None,
                script_file=None,
                ligand_params_files=(),
                spades_args=None,
                debug=False) :
    adopt_init_args(self, locals())
    assert (prefix is not None)
    self.xray_structure = fmodel.xray_structure.deep_copy_scatterers()
    self.r_work_start = fmodel.r_work()
    self.r_free_start = fmodel.r_free()
    self.best_r_work = 1.0
    self.best_r_free = 1.0
    self.worst_r_free = 0.0
    self.best_model = None
    self.best_hierarchy = None
    self.best_structure = None
    self.spades_args = spades_args
    self.d_min = fmodel.f_obs().d_min()
    fmodel.info().show_targets(out=out, text="starting model")
    self.data_file = prefix + "_data.mtz"
    self.pdb_file = prefix + "_in.pdb"
    fmodel.export_f_obs_flags_as_mtz(file_name=self.data_file)
    if (not os.path.isfile(self.pdb_file)) :
      open(self.pdb_file, "w").write(pdb_hierarchy.as_pdb_string(
        crystal_symmetry=fmodel.xray_structure))
    validation = validate(pdb_hierarchy=pdb_hierarchy,
      crystal_symmetry=fmodel.xray_structure, out=out,
      text="(starting model)")
    self.validation_start = validation
    self._process_info = get_process_info()
    self.setup_args()
    self.run_jobs()
    self.evaluate_models()

  def setup_args (self) :
    script_file = self.script_file
    if (script_file is None) :
      script_file = find_script_file("low_resolution_refine.xml")
      assert (script_file is not None)
    args = [
      "-parser:protocol %s" % script_file,
      "-s %s" % self.pdb_file, #os.path.abspath(self.pdb_file),
      "-mtzfile %s" % self.data_file, # os.path.abspath(self.data_file),
      "-run:preserve_header",
      "-crystal_refine",
      "-parser:script_vars symmdef=%s" % self.symm_file,
      "-parser:script_vars bfactstrat=%s" % self.params.adp_type,
      "-parser:script_vars map_type=%s" % self.params.map_type,
      "-nstruct 1",
    ]
    if (self.params.ignore_unrecognized) :
      args.append("-ignore_unrecognized_res")
    if (len(self.ligand_params_files) > 0) :
      args.append("-extra_res_fa %s" % " ".join(self.ligand_params_files))
    if (self.database is not None) :
      args.append("-database %s" % self.database)
    if (self.spades_args):
      with open(self.spades_args) as f:
        spades_args = f.read().splitlines()
      args += spades_args
    self.args = args
    if (self.debug) :
      print("", file=self.out)
      print("Rosetta command-line arguments:", file=self.out)
      for arg in self.args :
        print("  %s" % arg, file=self.out)
      print("", file=self.out)
    return args

  def run_jobs (self) :
    self.n_finished = 0
    self._process_info = get_process_info()
    nproc = self.mp_params.nproc
    all_args = []
    assert (self.params.number_of_models >= 1)
    n_jobs = self.params.number_of_models
    job_id = 1
    for i_model in range(n_jobs) :
      suffix = "_%d" % job_id
      pdb_out = self.prefix + "_in" + suffix + "_0001.pdb"
      job_args = [
        "-out:suffix %s" % suffix,
      ]
      args = list(self.args) + job_args
      all_args.append(rosetta_args(args=args, pdb_out=pdb_out,
        job_id=job_id))
      job_id += 1
    if (nproc in [None, Auto]) :
      if (self.mp_params.technology == "multiprocessing") :
        nproc = libtbx.introspection.number_of_processors()
      else :
        nproc = n_jobs
    nproc = min(nproc, n_jobs)
    print("Generating %d models on %d processors..." % \
      (n_jobs, nproc), file=self.out)
    print("", file=self.out)
    def check_result (result) :
      if (result.return_code > 0) :
        print("DEBUG dumping error message", file=self.out)
        print(result.error_message, file=self.out)
        raise RuntimeError("Aborted due to error in job %d." % result.job_id)
      elif (result.return_code < 0) :
        print("DEBUG dumping error message", file=self.out)
        print(result.error_message, file=self.out)
        raise RuntimeError("File '%s' not found." % result.file_name)
    if (nproc == 1) :
      _run = run_xtal_refine(self.debug)
      for _args in all_args :
        result = _run(_args)
        check_result(result)
        self.load_rosetta_model(result)
    else :
      from libtbx import easy_mp
      results = easy_mp.parallel_map(
        func=run_xtal_refine(debug=self.debug),
        iterable=all_args,
        params=self.mp_params)
      for result in results :
#        check_result(result)
        self.load_rosetta_model(result)

  def load_rosetta_model (self, result) :
    if (result.file_name is None) :
      print("ERROR in job %d" % result.job_id, file=self.out)
      print(result.error_message, file=self.out)
      return
    self.n_finished += 1
    from iotbx import file_reader
    new_pdb = file_reader.any_file(result.file_name, force_type="pdb")
    new_pdb.assert_file_type("pdb")
    new_hierarchy = new_pdb.file_object.hierarchy
    new_xrs = new_pdb.file_object.input.xray_structure_simple(
      crystal_symmetry=self.fmodel.xray_structure)
    pdb_atoms = self.pdb_hierarchy.atoms()
    score = None
    for line in new_pdb.file_object.input.remark_section() :
      if (line.startswith("REMARK   1 [ final ]")) :
        fields = line.strip().split()
        score = float(fields[-1].split("=")[-1])
        break
    try :
      rmsd = load_new_coordinates(
        old_atoms=pdb_atoms,
        new_atoms=new_hierarchy.atoms(),
        superpose=False,
        out=self.out)
    except Exception as e :
      print("Error calculating RMSD: %s" % str(e))
      rmsd = None
    self.fmodel.update_xray_structure(
      xray_structure=new_xrs,
      update_f_mask=True,
      update_f_calc=True)
    r_work = self.fmodel.r_work()
    r_free = self.fmodel.r_free()
    suffix = ""
    if (r_free < self.best_r_free) :
      self.best_r_free = r_free
      self.best_model = result.file_name
      self.best_hierarchy = new_hierarchy
      self.best_structure = new_xrs
      self.best_score = score
      self.best_rmsd = rmsd
      suffix = " ***"
    elif (r_free > self.worst_r_free) :
      self.worst_r_free = r_free
    fs = "  %3d: r_work = %6.4f  r_free = %6.4f  energy = %s  rmsd = %s%s"
    print(fs % (result.job_id, r_work, r_free,
      format_value("%8.2f", score), format_value("%5.3f", rmsd),  suffix), file=self.out)

  def evaluate_models (self) :
    self.fmodel.update_xray_structure(
      xray_structure=self.best_structure,
      update_f_mask=True,
      update_f_calc=True)
    print("", file=self.out)
    self.fmodel.info().show_targets(out=self.out, text="after ROSETTA")
    validate(self.best_hierarchy, self.fmodel.xray_structure, out=self.out, text="after ROSETTA")

def validate (pdb_hierarchy, crystal_symmetry, pdb_file=None, out=None, text="") :
  if (out is None) :
    out = sys.stdout
  from mmtbx.model import manager
  from mmtbx.validation import molprobity
  pdb_in = iotbx.pdb.input(
    source_info='string',
    lines=pdb_hierarchy.as_pdb_string(crystal_symmetry=crystal_symmetry))
  pdb_model = manager(pdb_in)
  validation = molprobity.molprobity(
    pdb_model,
    keep_hydrogens=False)
  if (text != "") : text = " " + text
  print("", file=out)
  print("  Validation statistics%s:" % text, file=out)
  validation.show_summary(out=out, prefix="    ")
  print("", file=out)
  return validation

class rosetta_args (object) :
  def __init__ (self, args, pdb_out, job_id) :
    adopt_init_args(self, locals())

class run_xtal_refine (object) :
  def __init__ (self, debug=False) :
    self.debug = debug
    self.env_phenix_rosetta_path = os.environ.get("PHENIX_ROSETTA_PATH", None)
    self.env_rosetta_bin = os.environ.get("ROSETTA_BIN", None)
    self.env_rosetta3_db = os.environ.get("ROSETTA3_DB", None)

  def __call__ (self, rosetta_args) :
    if ((self.env_phenix_rosetta_path is not None) and
        (not "PHENIX_ROSETTA_PATH" in os.environ)) :
      os.environ["PHENIX_ROSETTA_PATH"] = self.env_phenix_rosetta_path
    if ((self.env_rosetta_bin is not None) and
        (not "ROSETTA_BIN" in os.environ)) :
      os.environ["ROSETTA_BIN"] = self.env_rosetta_bin
    if ((self.env_rosetta3_db is not None) and
        (not "ROSETTA3_DB" in os.environ)) :
      os.environ["ROSETTA3_DB"] = self.env_rosetta3_db
    t1 = time.time()
    scripts_cmd = find_scripts_command()
    script_file = None # for debugging output
    if (self.debug):
      script_file = "rosetta_%s.sh" % rosetta_args.job_id
    result = run_hybrid_command(
      command_name=scripts_cmd,
      args=rosetta_args.args,
      out=null_out(),
      verbose=False,
      buffer_output=True,
      script_file=script_file)
    def save_debugging_info () :
      err_file = "job_%d_stderr.log" % rosetta_args.job_id
      out_file = "job_%d_stdout.log" % rosetta_args.job_id
      f1 = open(err_file, "w")
      f1.write("\n".join(result.stderr_lines))
      f2 = open(out_file, "w")
      f2.write("\n".join(result.stdout_lines))
    t2 = time.time()
    error_message = None
    pdb_out = rosetta_args.pdb_out
    return_code = result.return_code
    if (return_code != 0) :
      error_message = "ERROR: Rosetta exited with status %d\n" % \
        result.return_code
      error_message += "stderr output:\n" + "\n".join(result.stderr_lines)
      save_debugging_info()
      pdb_out = None
    elif (not os.path.isfile(pdb_out)) :
      return_code = -1
      error_message = "\n".join(result.stderr_lines)
      save_debugging_info()
    return rosetta_job_result(
      job_id=rosetta_args.job_id,
      command=scripts_cmd,
      process_time=(t2 - t1),
      file_name=pdb_out,
      error_message=error_message,
      return_code=return_code)

def run_hybrid_command (command_name, args, out=None, verbose=False,
    buffer_output=True, script_file=None):
  if (command_name is None) :
    raise RuntimeError("Couldn't find RosettaScripts binary on %s" %
      os.uname()[1])
  assert (command_name is not None)
  if (out is None) : out = sys.stdout
  build_path = abs(libtbx.env.build_path)
  phenix_path = libtbx.env.find_in_repositories(
    relative_path="phenix",
    test=os.path.isdir)
  assert (phenix_path is not None)
  # XXX avoid floating-point errors - is this necessary?
  #os.environ["BOOST_ADAPTBX_FPE_DEFAULT"] = "1"
  os.environ["ROSETTA_PHENIX_DIST"] = os.path.dirname(phenix_path)
  os.environ['LIBTBX_BUILD'] = os.path.join(
    os.path.dirname(os.environ['ROSETTA_PHENIX_DIST']),
    'build')
  os.environ["ROSETTA_PHENIX_BIN"] = sys.executable
  pythonpath = [ os.path.join(build_path, "lib") ]
  for module in libtbx.env.module_list :
    module_paths = module.assemble_pythonpath()
    for path in module_paths :
      if (not abs(path) in pythonpath) :
        pythonpath.append(abs(path))
  os.environ["ROSETTA_PHENIX_MODULES"] = ":".join(pythonpath)
  script_out = multi_out()
  if (script_file is not None) :
    f = open(script_file, "w")
    script_out.register("script", f)
  if (verbose) :
    script_out.register("stdout", sys.stdout)
  print("#" * 72, file=script_out)
  for env_name in ["ROSETTA_PHENIX_DIST", "ROSETTA_PHENIX_BIN",
                   "ROSETTA_PHENIX_MODULES", 'LIBTBX_BUILD'] :
    print("export %s=%s" % (env_name, os.environ[env_name]), file=script_out)
  if (sys.platform == "darwin") :
    dyld_name = "DYLD_FALLBACK_LIBRARY_PATH"
    dyld = os.environ.get(dyld_name)
    if dyld is None:
      print("export DYLD_LIBRARY_PATH=%s" % \
        os.environ["DYLD_LIBRARY_PATH"], file=script_out)
    else:
      print("export %s=%s" % (dyld_name, dyld), file=script_out)
  else :
    print("export LD_LIBRARY_PATH=%s" % \
      os.environ["LD_LIBRARY_PATH"], file=script_out)
  print("%s \\" % command_name, file=script_out)
  print("  %s" % " \\\n  ".join(args), file=script_out)
  print("#" * 72, file=script_out)
  if (not buffer_output):
    return easy_run.call(" ".join([command_name] + args))
  else :
    result = easy_run.fully_buffered(" ".join([command_name] + args))
    return result

def run_rosetta_scripts (args, debug=False) :
  pass

class rosetta_job_result (object) :
  def __init__ (self,
      job_id,
      command,
      process_time,
      return_code,
      file_name=None,
      error_message=None) :
    adopt_init_args(self, locals())

def load_new_coordinates (old_atoms, new_atoms, superpose=False, out=None) :
  if (out is not None) :
    out = sys.stdout
  old_atoms.reset_i_seq()
  new_atoms.reset_i_seq()
  from scitbx.array_family import flex
  atom_dict = {}
  # pass 1: select common atoms
  selection_old = flex.bool(old_atoms.size(), False)
  selection_new = flex.bool(new_atoms.size(), False)
  for i_seq, atom in enumerate(old_atoms) :
    atom_dict[atom.id_str()] = i_seq
  for j_seq, atom in enumerate(new_atoms) :
    id_str = atom.id_str()
    if (id_str in atom_dict) :
      i_seq = atom_dict[id_str]
      selection_old[i_seq] = True
      selection_new[j_seq] = True
  sites_all = new_atoms.extract_xyz()
  old_atoms_sel = old_atoms.select(selection_old)
  new_atoms_sel = new_atoms.select(selection_new)
  old_atoms_sel.reset_i_seq()
  new_atoms_sel.reset_i_seq()
  sites_sel = sites_all.select(selection_new)
  sites_old = old_atoms_sel.extract_xyz()
  # pass 2: determine permutation
  perm = flex.int(old_atoms_sel.size(), -1)
  for i_seq, atom in enumerate(old_atoms_sel) :
    atom_dict[atom.id_str()] = i_seq
  for j_seq, atom in enumerate(new_atoms_sel) :
    id_str = atom.id_str()
    if (id_str in atom_dict) :
      i_seq = atom_dict[id_str]
      perm[j_seq] = i_seq
  sites_sorted = flex.vec3_double(flex.select(sites_sel, permutation=perm))
  if (sites_sorted.size() != sites_old.size()) :
    raise RuntimeError("Mismatch in number of sites: %d vs. %d" %
      (sites_sorted.size(), sites_old.size()))
  # XXX this should be unnecessary if relaxing against electron density.  For
  # ordinary relax, the absolute positions may change arbitrarily, and will
  # need to be refit.  Probably should use SSM for this...
  if (superpose) :
    import scitbx.math.superpose
    lsq = scitbx.math.superpose.least_squares_fit(
      reference_sites=sites_old,
      other_sites=sites_sorted)
    sites_sorted = lsq.r.elems * sites_sorted + lsq.t.elems
    sites_all = lsq.r.elems * sites_all + lsq.t.elems
    new_atoms.set_xyz(flex.vec3_double(sites_all))
  rmsd = sites_old.rms_difference(sites_sorted)
  new_atoms.reset_i_seq()
  old_atoms.reset_i_seq()
  return rmsd

def rename_waters (pdb_hierarchy) :
  model = pdb_hierarchy.models()[0]
  for chain in model.chains() :
    for residue_group in chain.residue_groups() :
      for atom_group in residue_group.atom_groups() :
        if (atom_group.resname == "HOH") :
          atom_group.resname = "WAT"

def move_waters (pdb_hierarchy) :
  from iotbx.pdb import hierarchy
  model = pdb_hierarchy.models()[0]
  water_chains = []
  for chain in model.chains() :
    new_chain = hierarchy.chain()
    n_wats = 0
    new_chain.id = chain.id
    for residue_group in chain.residue_groups() :
      remove_rg = False
      for atom_group in residue_group.atom_groups() :
        if (atom_group.resname in ["HOH", "WAT"]) :
          atom_group_new = atom_group.detached_copy()
          residue_group.remove_atom_group(atom_group)
          residue_group_new = hierarchy.residue_group()
          residue_group_new.resseq = residue_group.resseq
          residue_group_new.icode = residue_group.icode
          residue_group_new.append_atom_group(atom_group_new)
          new_chain.append_residue_group(residue_group_new)
          n_wats += 1
          remove_rg = True
      if (remove_rg) :
        chain.remove_residue_group(residue_group)
    if (n_wats > 0) :
      water_chains.append(new_chain)
  if (len(water_chains) > 0) :
    for chain in water_chains :
      model.append_chain(chain)

def remove_heteroatoms (pdb_hierarchy, log) :
  if (log is None) : log = null_out()
  from iotbx.pdb import hierarchy
  assert (len(pdb_hierarchy.models()) == 1)
  model = pdb_hierarchy.models()[0]
  hetatm_chains = []
  for chain in model.chains() :
    n_hets = 0
    is_protein = chain.is_protein()
    new_chain = hierarchy.chain()
    for residue_group in chain.residue_groups() :
      for atom_group in residue_group.atom_groups() :
        resname = atom_group.resname
        if (resname in ["HOH", "WAT"]) :
          continue
        res_class = iotbx.pdb.common_residue_names_get_class(resname)
        if (not res_class in ["common_amino_acid", "common_element"]) :
          print("Warning: temporarily removing residue %s %s %s" % \
            (chain.id, residue_group.resid(), resname), file=log)
          atom_group_new = atom_group.detached_copy()
          residue_group.remove_atom_group(atom_group)
          residue_group_new = hierarchy.residue_group()
          residue_group_new.resseq = residue_group.resseq
          residue_group_new.icode = residue_group.icode
          residue_group_new.append_atom_group(atom_group_new)
          new_chain.append_residue_group(residue_group_new)
          n_hets += 1
      if (residue_group.atom_groups_size() == 0) :
        chain.remove_residue_group(residue_group)
    if (n_hets > 0) :
      hetatm_chains.append(new_chain)
    if (chain.residue_groups_size() == 0) :
      model.remove_chain(chain)
  return hetatm_chains
