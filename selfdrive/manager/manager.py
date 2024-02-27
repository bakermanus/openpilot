#!/usr/bin/env python3
import datetime
import os
import signal
import sys
import traceback
from typing import List, Tuple, Union

from cereal import log
import cereal.messaging as messaging
import openpilot.selfdrive.sentry as sentry
from openpilot.common.basedir import BASEDIR
from openpilot.common.params import Params, ParamKeyType
from openpilot.common.text_window import TextWindow
from openpilot.system.hardware import HARDWARE, PC
from openpilot.selfdrive.manager.helpers import unblock_stdout, write_onroad_params, save_bootlog
from openpilot.selfdrive.manager.process import ensure_running
from openpilot.selfdrive.manager.process_config import managed_processes
from openpilot.selfdrive.athena.registration import register, UNREGISTERED_DONGLE_ID
from openpilot.common.swaglog import cloudlog, add_file_handler
from openpilot.system.version import is_dirty, get_commit, get_version, get_origin, get_short_branch, \
                           get_normalized_origin, terms_version, training_version, \
                           is_tested_branch, is_release_branch, get_commit_date

from openpilot.selfdrive.frogpilot.functions.frogpilot_functions import DEFAULT_MODEL


def manager_init() -> None:
  save_bootlog()

  # Clear the error log on boot to prevent old errors from hanging around
  if os.path.isfile(os.path.join(sentry.CRASHES_DIR, 'error.txt')):
    os.remove(os.path.join(sentry.CRASHES_DIR, 'error.txt'))

  params = Params()
  params.clear_all(ParamKeyType.CLEAR_ON_MANAGER_START)
  params.clear_all(ParamKeyType.CLEAR_ON_ONROAD_TRANSITION)
  params.clear_all(ParamKeyType.CLEAR_ON_OFFROAD_TRANSITION)
  if is_release_branch():
    params.clear_all(ParamKeyType.DEVELOPMENT_ONLY)

  ############### Remove this after the April 26th update ###############

  previous_speed_limit = params.get_float("PreviousSpeedLimit")
  if previous_speed_limit >= 50:
    params.put_float("PreviousSpeedLimit", previous_speed_limit / 100)

  for priority_key in ["SLCPriority", "SLCPriority1", "SLCPriority2", "SLCPriority3"]:
    priority_value = params.get(priority_key)
    if isinstance(priority_value, int):
      params.remove(priority_key)

  attributes = ["AggressiveFollow", "StandardFollow", "RelaxedFollow", "AggressiveJerk", "StandardJerk", "RelaxedJerk"]
  values = {attr: params.get_float(attr) for attr in attributes}
  if any(value > 5 for value in values.values()):
    for attr, value in values.items():
      params.put_float(attr, value / 10)

  if params.get_bool("SilentMode"):
    attributes = ["DisengageVolume", "EngageVolume", "PromptVolume", "PromptDistractedVolume", "RefuseVolume", "WarningSoftVolume", "WarningImmediateVolume"]
    for attr in attributes:
      params.put_float(attr, 0)
    params.put_bool("SilentMode", False)

  #######################################################################

  # Check if the currently selected model still exists
  current_model = params.get("Model", encoding='utf-8')
  models_folder = os.path.join(BASEDIR, 'selfdrive/modeld/models/models')
  model_exists = current_model in [os.path.splitext(file)[0] for file in os.listdir(models_folder)]

  if not model_exists:
    params.remove("Model")

  FrogsGoMoo = HARDWARE.get_serial()[:3] == "cff"

  default_params: List[Tuple[str, Union[str, bytes]]] = [
    ("CompletedTrainingVersion", training_version if FrogsGoMoo else "0"),
    ("DisengageOnAccelerator", "0"),
    ("GsmMetered", "0" if FrogsGoMoo else "1"),
    ("HasAcceptedTerms", terms_version if FrogsGoMoo else "0"),
    ("LanguageSetting", "main_en"),
    ("OpenpilotEnabledToggle", "1"),
    ("LongitudinalPersonality", str(log.LongitudinalPersonality.standard)),

    # Default FrogPilot parameters
    ("AccelerationPath", "1"),
    ("AccelerationProfile", "3" if FrogsGoMoo else "2"),
    ("AdjacentPath", "1" if FrogsGoMoo else "0"),
    ("AdjacentPathMetrics", "1" if FrogsGoMoo else "0"),
    ("AdjustablePersonalities", "1"),
    ("AggressiveAcceleration", "1"),
    ("AggressiveFollow", "1" if FrogsGoMoo else "1.25"),
    ("AggressiveJerk", "0.55" if FrogsGoMoo else "0.5"),
    ("AlertVolumeControl", "0"),
    ("AlwaysOnLateral", "1"),
    ("AlwaysOnLateralMain", "1" if FrogsGoMoo else "0"),
    ("BlindSpotPath", "1"),
    ("CameraView", "1" if FrogsGoMoo else "0"),
    ("CECurves", "1"),
    ("CENavigation", "1"),
    ("CENavigationIntersections", "1"),
    ("CENavigationLead", "0" if FrogsGoMoo else "1"),
    ("CENavigationTurns", "1"),
    ("CESignal", "1"),
    ("CESlowerLead", "0"),
    ("CESpeed", "0"),
    ("CESpeedLead", "0"),
    ("CEStopLights", "1"),
    ("CEStopLightsLead", "0" if FrogsGoMoo else "1"),
    ("Compass", "1" if FrogsGoMoo else "0"),
    ("ConditionalExperimental", "1"),
    ("CrosstrekTorque", "0"),
    ("CurveSensitivity", "125" if FrogsGoMoo else "100"),
    ("CustomAlerts", "0"),
    ("CustomColors", "1"),
    ("CustomIcons", "1"),
    ("CustomPersonalities", "1"),
    ("CustomSignals", "1"),
    ("CustomSounds", "1"),
    ("CustomTheme", "1"),
    ("CustomUI", "1"),
    ("CydiaTune", "1"),
    ("DecelerationProfile", "1"),
    ("DeviceShutdown", "9"),
    ("DisableMTSCSmoothing", "0"),
    ("DisableVTSCSmoothing", "0"),
    ("DisengageVolume", "100"),
    ("DragonPilotTune", "0"),
    ("DriverCamera", "0"),
    ("DriveStats", "1"),
    ("DynamicPathWidth", "0"),
    ("EngageVolume", "100"),
    ("EVTable", "0" if FrogsGoMoo else "1"),
    ("ExperimentalModeActivation", "1"),
    ("ExperimentalModeViaLKAS", "1" if FrogsGoMoo else "0"),
    ("ExperimentalModeViaScreen", "0" if FrogsGoMoo else "1"),
    ("Fahrenheit", "0"),
    ("FireTheBabysitter", "1" if FrogsGoMoo else "0"),
    ("ForceAutoTune", "0"),
    ("FPSCounter", "1" if FrogsGoMoo else "0"),
    ("FrogsGoMooTune", "1" if FrogsGoMoo else "0"),
    ("FullMap", "0"),
    ("GasRegenCmd", "0"),
    ("GoatScream", "1"),
    ("GreenLightAlert", "0"),
    ("HideSpeed", "0"),
    ("HideSpeedUI", "0"),
    ("HigherBitrate", "1" if FrogsGoMoo else "0"),
    ("LaneChangeTime", "0"),
    ("LaneDetection", "1"),
    ("LaneDetectionWidth", "60"),
    ("LaneLinesWidth", "4"),
    ("LateralTune", "1"),
    ("LeadDepartingAlert", "0"),
    ("LeadInfo", "1" if FrogsGoMoo else "0"),
    ("LockDoors", "0"),
    ("LongitudinalTune", "1"),
    ("LongPitch", "0" if FrogsGoMoo else "1"),
    ("LoudBlindspotAlert", "0"),
    ("LowerVolt", "0" if FrogsGoMoo else "1"),
    ("MTSCAggressiveness", "100" if FrogsGoMoo else "100"),
    ("MTSCCurvatureCheck", "1" if FrogsGoMoo else "0"),
    ("MTSCLimit", "30" if FrogsGoMoo else "0"),
    ("Model", "certified-herbalist" if FrogsGoMoo else DEFAULT_MODEL),
    ("ModelUI", "1"),
    ("MTSCEnabled", "1"),
    ("MuteOverheated", "1" if FrogsGoMoo else "0"),
    ("NNFF", "1"),
    ("NoLogging", "0"),
    ("NoUploads", "0"),
    ("NudgelessLaneChange", "1"),
    ("NumericalTemp", "1" if FrogsGoMoo else "0"),
    ("Offset1", "5"),
    ("Offset2", "7" if FrogsGoMoo else "5"),
    ("Offset3", "10" if FrogsGoMoo else "5"),
    ("Offset4", "20" if FrogsGoMoo else "10"),
    ("OneLaneChange", "1"),
    ("PathEdgeWidth", "20"),
    ("PathWidth", "61"),
    ("PauseLateralOnSignal", "0"),
    ("PedalsOnUI", "1"),
    ("PersonalitiesViaScreen", "0" if FrogsGoMoo else "1"),
    ("PersonalitiesViaWheel", "1"),
    ("PreferredSchedule", "0"),
    ("PromptVolume", "100"),
    ("PromptDistractedVolume", "100"),
    ("QOLControls", "1"),
    ("QOLVisuals", "1"),
    ("RandomEvents", "1" if FrogsGoMoo else "0"),
    ("RefuseVolume", "100"),
    ("RelaxedFollow", "3.0" if FrogsGoMoo else "1.75"),
    ("RelaxedJerk", "5.0" if FrogsGoMoo else "1.0"),
    ("ReverseCruise", "1" if FrogsGoMoo else "0"),
    ("ReverseCruiseUI", "1"),
    ("RoadEdgesWidth", "2"),
    ("RoadNameUI", "1"),
    ("RotatingWheel", "1"),
    ("ScreenBrightness", "101"),
    ("SearchInput", "0"),
    ("SetSpeedLimit", "0"),
    ("SetSpeedOffset", "0"),
    ("ShowCPU", "1" if FrogsGoMoo else "0"),
    ("ShowGPU", "0"),
    ("ShowIP", "0"),
    ("ShowMemoryUsage", "1" if FrogsGoMoo else "0"),
    ("Sidebar", "1" if FrogsGoMoo else "0"),
    ("SLCConfirmation", "1"),
    ("SLCConfirmationLower", "1"),
    ("SLCConfirmationHigher", "0"),
    ("SLCFallback", "2"),
    ("SLCOverride", "1"),
    ("SLCPriority1", "Dashboard"),
    ("SLCPriority2", "Offline Maps"),
    ("SLCPriority3", "Navigation"),
    ("SmoothBraking", "1"),
    ("SNGHack", "0" if FrogsGoMoo else "1"),
    ("SpeedLimitChangedAlert", "0"),
    ("SpeedLimitController", "1"),
    ("StandardFollow", "1.45"),
    ("StandardJerk", "1.0"),
    ("StoppingDistance", "3" if FrogsGoMoo else "0"),
    ("TurnAggressiveness", "150" if FrogsGoMoo else "100"),
    ("TurnDesires", "1" if FrogsGoMoo else "0"),
    ("UnlimitedLength", "1"),
    ("UseLateralJerk", "0"),
    ("UseSI", "1" if FrogsGoMoo else "0"),
    ("UseVienna", "0"),
    ("VisionTurnControl", "1"),
    ("WarningSoftVolume", "100"),
    ("WarningImmediateVolume", "100"),
    ("WheelIcon", "1" if FrogsGoMoo else "3"),
    ("WheelSpeed", "0")
  ]
  if not PC:
    default_params.append(("LastUpdateTime", datetime.datetime.utcnow().isoformat().encode('utf8')))

  if params.get_bool("RecordFrontLock"):
    params.put_bool("RecordFront", True)

  # set unset params
  retry_params = []
  for k, v in default_params:
    if params.get(k) is None:
      try:
        params.put(k, v)
      except Exception:
        retry_params.append((k, v))

  for k, v in retry_params:
    try:
      if params.get(k) is None:
        params.put(k, v)
    except Exception:
      pass

  # Create folders needed for msgq
  try:
    os.mkdir("/dev/shm")
  except FileExistsError:
    pass
  except PermissionError:
    print("WARNING: failed to make /dev/shm")

  # set version params
  params.put("Version", get_version())
  params.put("TermsVersion", terms_version)
  params.put("TrainingVersion", training_version)
  params.put("GitCommit", get_commit())
  params.put("GitCommitDate", get_commit_date())
  params.put("GitBranch", get_short_branch())
  params.put("GitRemote", get_origin())
  params.put_bool("IsTestedBranch", is_tested_branch())
  params.put_bool("IsReleaseBranch", is_release_branch())

  # set dongle id
  reg_res = register(show_spinner=True)
  if reg_res:
    dongle_id = reg_res
  else:
    serial = params.get("HardwareSerial")
    raise Exception(f"Registration failed for device {serial}")
  os.environ['DONGLE_ID'] = dongle_id  # Needed for swaglog
  os.environ['GIT_ORIGIN'] = get_normalized_origin() # Needed for swaglog
  os.environ['GIT_BRANCH'] = get_short_branch() # Needed for swaglog
  os.environ['GIT_COMMIT'] = get_commit() # Needed for swaglog

  if not is_dirty():
    os.environ['CLEAN'] = '1'

  # init logging
  sentry.init(sentry.SentryProject.SELFDRIVE)
  cloudlog.bind_global(dongle_id=dongle_id,
                       version=get_version(),
                       origin=get_normalized_origin(),
                       branch=get_short_branch(),
                       commit=get_commit(),
                       dirty=is_dirty(),
                       device=HARDWARE.get_device_type())

  # preimport all processes
  for p in managed_processes.values():
    p.prepare()


def manager_cleanup() -> None:
  # send signals to kill all procs
  for p in managed_processes.values():
    p.stop(block=False)

  # ensure all are killed
  for p in managed_processes.values():
    p.stop(block=True)

  cloudlog.info("everything is dead")


def update_frogpilot_params(params, params_memory):
  keys = ["DisableOnroadUploads", "FireTheBabysitter", "NoLogging", "NoUploads", "RoadNameUI"]
  for key in keys:
    params_memory.put_bool(key, params.get_bool(key))

def manager_thread() -> None:
  cloudlog.bind(daemon="manager")
  cloudlog.info("manager start")
  cloudlog.info({"environ": os.environ})

  params = Params()
  params_memory = Params("/dev/shm/params")

  update_frogpilot_params(params, params_memory)

  ignore: List[str] = []
  if params.get("DongleId", encoding='utf8') in (None, UNREGISTERED_DONGLE_ID):
    ignore += ["manage_athenad", "uploader"]
  if os.getenv("NOBOARD") is not None:
    ignore.append("pandad")
  ignore += [x for x in os.getenv("BLOCK", "").split(",") if len(x) > 0]

  sm = messaging.SubMaster(['deviceState', 'carParams'], poll='deviceState')
  pm = messaging.PubMaster(['managerState'])

  write_onroad_params(False, params)
  ensure_running(managed_processes.values(), False, params=params, params_memory=params_memory, CP=sm['carParams'], not_run=ignore)

  started_prev = False

  while True:
    sm.update(1000)

    started = sm['deviceState'].started

    if started and not started_prev:
      params.clear_all(ParamKeyType.CLEAR_ON_ONROAD_TRANSITION)
    elif not started and started_prev:
      params.clear_all(ParamKeyType.CLEAR_ON_OFFROAD_TRANSITION)

      # Clear the error log on offroad transition to prevent old errors from hanging around
      if os.path.isfile(os.path.join(sentry.CRASHES_DIR, 'error.txt')):
        os.remove(os.path.join(sentry.CRASHES_DIR, 'error.txt'))

    # update onroad params, which drives boardd's safety setter thread
    if started != started_prev:
      write_onroad_params(started, params)

    started_prev = started

    ensure_running(managed_processes.values(), started, params=params, params_memory=params_memory, CP=sm['carParams'], not_run=ignore)

    running = ' '.join("%s%s\u001b[0m" % ("\u001b[32m" if p.proc.is_alive() else "\u001b[31m", p.name)
                       for p in managed_processes.values() if p.proc)
    print(running)
    cloudlog.debug(running)

    # send managerState
    msg = messaging.new_message('managerState', valid=True)
    msg.managerState.processes = [p.get_process_state_msg() for p in managed_processes.values()]
    pm.send('managerState', msg)

    # Exit main loop when uninstall/shutdown/reboot is needed
    shutdown = False
    for param in ("DoUninstall", "DoShutdown", "DoReboot", "DoSoftReboot"):
      if params.get_bool(param):
        shutdown = True
        params.put("LastManagerExitReason", f"{param} {datetime.datetime.now()}")
        cloudlog.warning(f"Shutting down manager - {param} set")

    if shutdown:
      break

    if params_memory.get_bool("FrogPilotTogglesUpdated"):
      update_frogpilot_params(params, params_memory)

def main() -> None:
  manager_init()
  if os.getenv("PREPAREONLY") is not None:
    return

  # SystemExit on sigterm
  signal.signal(signal.SIGTERM, lambda signum, frame: sys.exit(1))

  try:
    manager_thread()
  except Exception:
    traceback.print_exc()
    sentry.capture_exception()
  finally:
    manager_cleanup()

  params = Params()
  if params.get_bool("DoUninstall"):
    cloudlog.warning("uninstalling")
    HARDWARE.uninstall()
  elif params.get_bool("DoSoftReboot"):
    cloudlog.warning("softreboot")
    HARDWARE.soft_reboot()
  elif params.get_bool("DoReboot"):
    cloudlog.warning("reboot")
    HARDWARE.reboot()
  elif params.get_bool("DoShutdown"):
    cloudlog.warning("shutdown")
    HARDWARE.shutdown()


if __name__ == "__main__":
  unblock_stdout()

  try:
    main()
  except KeyboardInterrupt:
    print("got CTRL-C, exiting")
  except Exception:
    add_file_handler(cloudlog)
    cloudlog.exception("Manager failed to start")

    try:
      managed_processes['ui'].stop()
    except Exception:
      pass

    # Show last 3 lines of traceback
    error = traceback.format_exc(-3)
    error = "Manager failed to start\n\n" + error
    with TextWindow(error) as t:
      t.wait_for_exit()

    raise

  # manual exit because we are forked
  sys.exit(0)
