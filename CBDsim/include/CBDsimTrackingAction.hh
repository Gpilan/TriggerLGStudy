#ifndef CBDsimTrackingAction_h
#define CBDsimTrackingAction_h 1

#include "G4UserTrackingAction.hh"

class CBDsimTrackingAction : public G4UserTrackingAction {
public:
  CBDsimTrackingAction();
  virtual ~CBDsimTrackingAction();

  virtual void PostUserTrackingAction(const G4Track* track) override;
};

#endif
