// Included only by the opt-in CBDsimPrecision target, inside diagnostics namespace.
// Read-only observer: no random draws, track mutations, or per-step flushes.
struct PrecisionTrack {
  G4ThreeVector birth;
  G4double birthTime = 0., birthEnergy = 0., length = 0.;
  G4String birthPath, lastPrePath = "none", lastPostPath = "none";
  G4double lastTime = 0., lastEnergy = 0.;
  G4ThreeVector lastDirection;
  G4long steps = 0, toWorld = 0, fromWorld = 0;
  G4long shoulderReflections = 0, shoulderAbsorptions = 0;
  G4String firstWorld = "none", lastWorld = "none";
};

struct PrecisionTrace {
  bool enabled = std::getenv("CBDsim_PRECISION_OUT") != nullptr;
  G4long maxTracks = Limit("CBDsim_PRECISION_MAX_TRACKS", 100000);
  G4long maxSteps = Limit("CBDsim_PRECISION_MAX_STEPS", 0);
  G4long eventFilter = Limit("CBDsim_PRECISION_EVENT", -1);
  G4long trackFilter = Limit("CBDsim_PRECISION_TRACK", -1);
  G4long tracksRecorded = 0, stepsRecorded = 0, stepsOmitted = 0;
  G4long eventSeen = 0, eventSelected = 0, eventOmitted = 0;
  std::unordered_map<G4int, PrecisionTrack> tracks;
  std::ofstream photons, steps, events;

  static G4long Limit(const char* name, G4long fallback) {
    const char* v = std::getenv(name);
    if (!v) return fallback;
    char* end = nullptr;
    const auto n = std::strtol(v, &end, 10);
    if (!*v || *end || n < -1) {
      G4Exception("PrecisionTrace", "PrecisionConfig", FatalException, name);
    }
    return n;
  }
  static G4String Path(const G4VTouchable* touch) {
    if (!touch || !touch->GetVolume()) return "none";
    std::ostringstream out;
    for (int i = touch->GetHistoryDepth(); i >= 0; --i)
      out << '/' << touch->GetVolume(i)->GetName() << '[' << touch->GetCopyNumber(i) << ']';
    return out.str();
  }
  static void Vec(std::ostream& out, const G4ThreeVector& p, double unit = CLHEP::mm) {
    out << p.x()/unit << '\t' << p.y()/unit << '\t' << p.z()/unit;
  }
  static G4String Crossing(const G4Step* step) {
    // Compact field: step|source touchable|world xyz|source-local xyz (mm).
    const auto* pre = step->GetPreStepPoint();
    const auto p = step->GetPostStepPoint()->GetPosition();
    auto local = p;
    if (const auto* touch = dynamic_cast<const G4TouchableHistory*>(pre->GetTouchable()))
      local = touch->GetHistory()->GetTopTransform().TransformPoint(p);
    std::ostringstream out;
    out << std::setprecision(17) << step->GetTrack()->GetCurrentStepNumber() << '|'
        << Path(pre->GetTouchable()) << '|' << p.x()/CLHEP::mm << ',' << p.y()/CLHEP::mm
        << ',' << p.z()/CLHEP::mm << '|' << local.x()/CLHEP::mm << ','
        << local.y()/CLHEP::mm << ',' << local.z()/CLHEP::mm;
    return out.str();
  }
  void Open() {
    if (photons.is_open()) return;
    const auto base = G4String(std::getenv("CBDsim_PRECISION_OUT")) + "/";
    const auto suffix = "-w" + std::to_string(G4Threading::G4GetThreadId()) + ".tsv";
    // The runner creates a fresh directory. Stream failures are fatal, never silent.
    photons.exceptions(std::ios::badbit | std::ios::failbit);
    steps.exceptions(std::ios::badbit | std::ios::failbit);
    events.exceptions(std::ios::badbit | std::ios::failbit);
    photons.open(base + "photons" + suffix);
    steps.open(base + "steps" + suffix);
    events.open(base + "events" + suffix);
    photons << std::setprecision(17);
    steps << std::setprecision(17);
    photons << "event\ttrack\tparent\tsource\tbirth_path\tbirth_x_mm\tbirth_y_mm\tbirth_z_mm"
      "\tbirth_time_ns\tbirth_energy_eV\tfate\tfinal_region\tfinal_process\ttrack_status"
      "\tfinal_pre_volume\tfinal_post_volume\tfinal_material\tfinal_boundary_status"
      "\tfinal_x_mm\tfinal_y_mm\tfinal_z_mm\tpath_mm\tsteps\tto_world_transmissions"
      "\tfrom_world_transmissions\tfirst_to_world\tlast_to_world\tevent_aborted"
      "\tlast_pre_path\tlast_post_path\tlast_time_ns\tlast_energy_eV\tlast_dx\tlast_dy\tlast_dz\tshoulder_reflections\tshoulder_absorptions\n";
    steps << "event\ttrack\tstep\tpre_path\tpost_path\tpre_material\tpost_material"
      "\tpre_x_mm\tpre_y_mm\tpre_z_mm\tpost_x_mm\tpost_y_mm\tpost_z_mm"
      "\tpre_dx\tpre_dy\tpre_dz\tpost_dx\tpost_dy\tpost_dz"
      "\tpre_time_ns\tpost_time_ns\tpre_energy_eV\tpost_energy_eV\tlength_mm"
      "\tprocess\tstep_status\tboundary_status\ttrack_status\n";
    events << "event\taborted\tseen\tselected\tomitted_by_track_cap\ttracks_written"
      "\tsteps_written_total\tsteps_omitted_by_cap_total\n";
  }
  void BeginTrack(const G4Track* track) {
    if (!enabled) return;
    Open();
    ++eventSeen;
    if (eventFilter >= 0 && Budget().id != eventFilter) return;
    ++eventSelected;
    if (tracksRecorded >= maxTracks) { ++eventOmitted; return; }
    ++tracksRecorded;
    auto& t = tracks[track->GetTrackID()];
    t.birth = track->GetPosition();
    t.birthTime = track->GetGlobalTime();
    t.birthEnergy = track->GetKineticEnergy();
    t.birthPath = Path(track->GetTouchable());
  }
  void Step(const G4Step* step) {
    if (!enabled) return;
    const auto* track = step->GetTrack();
    auto it = tracks.find(track->GetTrackID());
    if (it == tracks.end()) return;
    auto& t = it->second;
    ++t.steps;
    t.length += step->GetStepLength();
    const auto* pre = step->GetPreStepPoint();
    const auto* post = step->GetPostStepPoint();
    if (track->GetTrackStatus() == fStopAndKill || track->GetTrackStatus() == fKillTrackAndSecondaries) {
      t.lastPrePath = Path(pre->GetTouchable());
      t.lastPostPath = Path(post->GetTouchable());
      t.lastTime = post->GetGlobalTime();
      t.lastEnergy = post->GetKineticEnergy();
      t.lastDirection = post->GetMomentumDirection();
    }
    int boundary = -1;
    if (post->GetStepStatus() == fGeomBoundary)
      if (auto* b = FindOpBoundaryProcess(track)) boundary = static_cast<int>(b->GetStatus());
    // Diagnostic for the current 0.2 mm round-tip model: only internal LG encounters
    // with its exposed cylinder-start plane, not zero-length relocation steps.
    if ((boundary == SpikeReflection || boundary == Absorption) && pre->GetPhysicalVolume()
        && post->GetPhysicalVolume() && pre->GetPhysicalVolume()->GetName()=="protoLightGuidePhys"
        && post->GetPhysicalVolume()->GetName()=="protoWorldPhys") {
      if (const auto* touch = dynamic_cast<const G4TouchableHistory*>(pre->GetTouchable())) {
        const auto local = touch->GetHistory()->GetTopTransform().TransformPoint(post->GetPosition());
        if (std::abs(local.y()/CLHEP::mm + 59.82) < 1.e-6
            && std::hypot(local.x(),local.z())/CLHEP::mm <= 7.500001) {
          if (boundary == SpikeReflection) ++t.shoulderReflections;
          else ++t.shoulderAbsorptions;
        }
      }
    }
    if (boundary == Transmission || boundary == FresnelRefraction || boundary == SameMaterial) {
      const auto a = VolRegion(pre->GetPhysicalVolume() ? pre->GetPhysicalVolume()->GetName() : "none");
      const auto b = VolRegion(post->GetPhysicalVolume() ? post->GetPhysicalVolume()->GetName() : "none");
      if (a != "world" && b == "world") {
        ++t.toWorld;
        t.lastWorld = Crossing(step);
        if (t.toWorld == 1) t.firstWorld = t.lastWorld;
      }
      if (a == "world" && b != "world") ++t.fromWorld;
    }
    if (maxSteps == 0 || (trackFilter >= 0 && track->GetTrackID() != trackFilter)) return;
    if (stepsRecorded >= maxSteps) { ++stepsOmitted; return; }
    ++stepsRecorded;
    steps << Budget().id << '\t' << track->GetTrackID() << '\t' << track->GetCurrentStepNumber()
      << '\t' << Path(pre->GetTouchable()) << '\t' << Path(post->GetTouchable())
      << '\t' << (pre->GetMaterial() ? pre->GetMaterial()->GetName() : "none")
      << '\t' << (post->GetMaterial() ? post->GetMaterial()->GetName() : "none") << '\t';
    Vec(steps, pre->GetPosition()); steps << '\t'; Vec(steps, post->GetPosition()); steps << '\t';
    Vec(steps, pre->GetMomentumDirection(), 1.); steps << '\t';
    Vec(steps, post->GetMomentumDirection(), 1.);
    steps << '\t' << pre->GetGlobalTime()/CLHEP::ns << '\t' << post->GetGlobalTime()/CLHEP::ns
      << '\t' << pre->GetKineticEnergy()/CLHEP::eV << '\t' << post->GetKineticEnergy()/CLHEP::eV
      << '\t' << step->GetStepLength()/CLHEP::mm << '\t'
      << (post->GetProcessDefinedStep() ? post->GetProcessDefinedStep()->GetProcessName() : "none")
      << '\t' << static_cast<int>(post->GetStepStatus()) << '\t' << boundary
      << '\t' << static_cast<int>(track->GetTrackStatus()) << '\n';
  }
  void EndEvent(bool aborted) {
    if (!enabled) return;
    Open();
    for (const auto& entry : tracks) {
      const auto& t = entry.second;
      const auto& b = Budget().tracks.at(entry.first);
      photons << Budget().id << '\t' << entry.first << '\t' << b.parent << '\t' << b.source
        << '\t' << t.birthPath << '\t';
      Vec(photons, t.birth);
      photons << '\t' << t.birthTime/CLHEP::ns << '\t' << t.birthEnergy/CLHEP::eV
        << '\t' << (b.terminal.empty() ? "unfinished" : b.terminal) << '\t' << b.region
        << '\t' << b.process << '\t' << b.status << '\t' << b.preVolume << '\t' << b.postVolume
        << '\t' << b.postMaterial << '\t' << b.boundaryStatus << '\t';
      Vec(photons, b.position);
      photons << '\t' << t.length/CLHEP::mm << '\t' << t.steps << '\t' << t.toWorld
        << '\t' << t.fromWorld << '\t' << t.firstWorld << '\t' << t.lastWorld << '\t' << aborted
        << '\t' << t.lastPrePath << '\t' << t.lastPostPath << '\t' << t.lastTime/CLHEP::ns
        << '\t' << t.lastEnergy/CLHEP::eV << '\t';
      Vec(photons, t.lastDirection, 1.); photons << '\t' << t.shoulderReflections << '\t' << t.shoulderAbsorptions << '\n';
    }
    events << Budget().id << '\t' << aborted << '\t' << eventSeen << '\t' << eventSelected
      << '\t' << eventOmitted << '\t' << tracks.size() << '\t' << stepsRecorded << '\t' << stepsOmitted << '\n';
    photons.flush(); steps.flush(); events.flush();
    tracks.clear(); eventSeen = eventSelected = eventOmitted = 0;
  }
};
thread_local std::unique_ptr<PrecisionTrace> precision;
PrecisionTrace& Precision() {
  if (!precision) precision = std::make_unique<PrecisionTrace>();
  return *precision;
}
void ResetPrecision() { precision.reset(); }
