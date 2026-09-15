#include "CBDsimTrackingAction.hh"
#include "CBDsimOpticalDiagnostics.hh"

CBDsimTrackingAction::CBDsimTrackingAction() : G4UserTrackingAction() {}

CBDsimTrackingAction::~CBDsimTrackingAction() {}

void CBDsimTrackingAction::PostUserTrackingAction(const G4Track* track) {
  CBDsimOpticalDiagnostics::PostUserTrackingAction(track);
}

void CBDsimTrackingAction::PreUserTrackingAction(const G4Track* track) {
  CBDsimOpticalDiagnostics::PreUserTrackingAction(track);
}
