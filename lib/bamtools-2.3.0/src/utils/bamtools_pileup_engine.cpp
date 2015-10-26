// ***************************************************************************
// bamtools_pileup_engine.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 9 March 2012 (DB)
// ---------------------------------------------------------------------------
// Provides pileup at position functionality for various tools.
// ***************************************************************************

#include "utils/bamtools_pileup_engine.h"
using namespace BamTools;

#include <iostream>
using namespace std;

// ---------------------------------------------
// PileupEnginePrivate implementation

struct PileupEngine::PileupEnginePrivate {
  
    // data members
    int CurrentId;
    int CurrentPosition;
    vector<BamAlignment> CurrentAlignments;
    PileupPosition CurrentPileupData;
    
    bool IsFirstAlignment;
    vector<PileupVisitor*> Visitors;
  
    // ctor & dtor
    PileupEnginePrivate(void)
        : CurrentId(-1)
        , CurrentPosition(-1)
        , IsFirstAlignment(true)
    { }
    ~PileupEnginePrivate(void) { }
    
    // 'public' methods
    bool AddAlignment(const BamAlignment& al);
    void Flush(void);
    
    // internal methods
    private:
        void ApplyVisitors(void);
        void ClearOldData(void);
        void CreatePileupData(void);
        void ParseAlignmentCigar(const BamAlignment& al);
};

bool PileupEngine::PileupEnginePrivate::AddAlignment(const BamAlignment& al) {
  
    // if first time
    if ( IsFirstAlignment ) {
      
        // set initial markers 
        CurrentId       = al.RefID;
        CurrentPosition = al.Position;
        
        // store first entry
        CurrentAlignments.clear();
        CurrentAlignments.push_back(al);
        
        // set flag & return
        IsFirstAlignment = false;
        return true;
    }
  
    // if same reference
    if ( al.RefID == CurrentId ) {
      
        // if same position, store and move on
        if ( al.Position == CurrentPosition )
            CurrentAlignments.push_back(al);
        
        // if less than CurrentPosition - sorting error => ABORT
        else if ( al.Position < CurrentPosition ) {
            cerr << "Pileup::Run() : Data not sorted correctly!" << endl;
            return false;
        }
        
        // else print pileup data until 'catching up' to CurrentPosition
        else {
            while ( al.Position > CurrentPosition ) {
                ApplyVisitors();
                ++CurrentPosition;
            }
            CurrentAlignments.push_back(al);
        }
    } 

    // if reference ID less than CurrentId - sorting error => ABORT
    else if ( al.RefID < CurrentId ) {
        cerr << "Pileup::Run() : Data not sorted correctly!" << endl;
        return false;
    }

    // else moved forward onto next reference
    else {
        
        // print any remaining pileup data from previous reference
        while ( !CurrentAlignments.empty() ) {
            ApplyVisitors();
            ++CurrentPosition;
        }
        
        // store first entry on this new reference, update markers
        CurrentAlignments.clear();
        CurrentAlignments.push_back(al);
        CurrentId = al.RefID;
        CurrentPosition = al.Position;
    }
  
    return true;
}

void PileupEngine::PileupEnginePrivate::ApplyVisitors(void) {
  
    // parse CIGAR data in BamAlignments to build up current pileup data
    CreatePileupData();
  
    // apply all visitors to current alignment set
    vector<PileupVisitor*>::const_iterator visitorIter = Visitors.begin();
    vector<PileupVisitor*>::const_iterator visitorEnd  = Visitors.end();
    for ( ; visitorIter != visitorEnd; ++visitorIter ) 
        (*visitorIter)->Visit(CurrentPileupData);
}

void PileupEngine::PileupEnginePrivate::ClearOldData(void) {
 
    // remove any alignments that end before our CurrentPosition
    // N.B. - BAM positions are 0-based, half-open. GetEndPosition() returns a 1-based position,
    //        while our CurrentPosition is 0-based. For example, an alignment with 'endPosition' of
    //        100 does not overlap a 'CurrentPosition' of 100, and should be discarded.

    size_t i = 0;
    size_t j = 0;
    const size_t numAlignments = CurrentAlignments.size();
    while ( i < numAlignments ) {

        // skip over alignment if its (1-based) endPosition is <= to (0-based) CurrentPosition
        // i.e. this entry will not be saved upon vector resize
        const int endPosition = CurrentAlignments[i].GetEndPosition();
        if ( endPosition <= CurrentPosition ) {
            ++i;
            continue;
        }

        // otherwise alignment ends after CurrentPosition
        // move it towards vector beginning, at index j
        if ( i != j )
            CurrentAlignments[j] = CurrentAlignments[i];

        // increment our indices
        ++i;
        ++j;
    }

    // 'squeeze' vector to size j, discarding all remaining alignments in the container
    CurrentAlignments.resize(j);
}

void PileupEngine::PileupEnginePrivate::CreatePileupData(void) {
  
    // remove any non-overlapping alignments
    ClearOldData();
  
    // set pileup refId, position to current markers
    CurrentPileupData.RefId = CurrentId;
    CurrentPileupData.Position = CurrentPosition;
    CurrentPileupData.PileupAlignments.clear();
    
    // parse CIGAR data in remaining alignments 
    vector<BamAlignment>::const_iterator alIter = CurrentAlignments.begin();
    vector<BamAlignment>::const_iterator alEnd  = CurrentAlignments.end(); 
    for ( ; alIter != alEnd; ++alIter )
        ParseAlignmentCigar( (*alIter) );
}

void PileupEngine::PileupEnginePrivate::Flush(void) {
    while ( !CurrentAlignments.empty() ) {
        ApplyVisitors();
        ++CurrentPosition;
    }
}

void PileupEngine::PileupEnginePrivate::ParseAlignmentCigar(const BamAlignment& al) {
  
    // skip if unmapped
    if ( !al.IsMapped() ) return;
    
    // intialize local variables
    int  genomePosition      = al.Position;
    int  positionInAlignment = 0;
    bool isNewReadSegment    = true;
    bool saveAlignment       = true;    
    PileupAlignment pileupAlignment(al);
    
    // iterate over CIGAR operations
    const int numCigarOps = (const int)al.CigarData.size();
    for (int i = 0; i < numCigarOps; ++i ) { 
        const CigarOp& op = al.CigarData.at(i);
      
        // if op is MATCH
        if ( op.Type == 'M' ) {
          
            // if match op overlaps current position
            if ( genomePosition + (int)op.Length > CurrentPosition ) {
              
                // set pileup data
                pileupAlignment.IsCurrentDeletion   = false;
                pileupAlignment.IsNextDeletion      = false;
                pileupAlignment.IsNextInsertion     = false;
                pileupAlignment.PositionInAlignment = positionInAlignment + (CurrentPosition - genomePosition);
                
                // check for beginning of read segment
                if ( genomePosition == CurrentPosition && isNewReadSegment ) 
                    pileupAlignment.IsSegmentBegin = true;
                
                // if we're at the end of a match operation
                if ( genomePosition + (int)op.Length - 1 == CurrentPosition ) {
                    
                    // if not last operation
                    if ( i < numCigarOps - 1 ) {
                        
                        // check next CIGAR op
                        const CigarOp& nextOp = al.CigarData.at(i+1);
                        
                        // if next CIGAR op is DELETION
                        if ( nextOp.Type == 'D') {
                            pileupAlignment.IsNextDeletion = true;
                            pileupAlignment.DeletionLength = nextOp.Length;
                        }
                        
                        // if next CIGAR op is INSERTION
                        else if ( nextOp.Type == 'I' ) {
                            pileupAlignment.IsNextInsertion = true;
                            pileupAlignment.InsertionLength = nextOp.Length;
                        }
                            
                        // if next CIGAR op is either DELETION or INSERTION
                        if ( nextOp.Type == 'D' || nextOp.Type == 'I' ) {

                            // if there is a CIGAR op after the DEL/INS
                            if ( i < numCigarOps - 2 ) {
                                const CigarOp& nextNextOp = al.CigarData.at(i+2);
                                
                                // if next CIGAR op is clipping or ref_skip
                                if ( nextNextOp.Type == 'S' || 
                                     nextNextOp.Type == 'N' ||
                                     nextNextOp.Type == 'H' )
                                    pileupAlignment.IsSegmentEnd = true;
                            } 
                            else {
                                pileupAlignment.IsSegmentEnd = true;
                                
                                // if next CIGAR op is clipping or ref_skip
                                if ( nextOp.Type == 'S' || 
                                     nextOp.Type == 'N' ||
                                     nextOp.Type == 'H' )
                                    pileupAlignment.IsSegmentEnd = true;
                            }
                        }
                        
                        // otherwise
                        else { 
                        
                            // if next CIGAR op is clipping or ref_skip
                            if ( nextOp.Type == 'S' || 
                                 nextOp.Type == 'N' ||
                                 nextOp.Type == 'H' )
                                pileupAlignment.IsSegmentEnd = true;
                        }
                    }
                    
                    // else this is last operation
                    else pileupAlignment.IsSegmentEnd = true;
                }
            }
          
            // increment markers
            genomePosition      += op.Length;
            positionInAlignment += op.Length;
        } 
        
        // if op is DELETION
        else if ( op.Type == 'D' ) {
          
            // if deletion op overlaps current position
            if ( genomePosition + (int)op.Length > CurrentPosition ) {
              
                // set pileup data
                pileupAlignment.IsCurrentDeletion   = true;
                pileupAlignment.IsNextDeletion      = false;
                pileupAlignment.IsNextInsertion     = true;
                pileupAlignment.PositionInAlignment = positionInAlignment + (CurrentPosition - genomePosition);
            }
            
            // increment marker
            genomePosition += op.Length;
        }

        // if op is REF_SKIP
        else if ( op.Type == 'N' ) {
            genomePosition += op.Length;
        }
        
        // if op is INSERTION or SOFT_CLIP
        else if ( op.Type == 'I' || op.Type == 'S' ) {
            positionInAlignment += op.Length;
        }
        
        // checl for beginning of new read segment
        if ( op.Type == 'N' ||
             op.Type == 'S' ||
             op.Type == 'H' )
            isNewReadSegment = true;
        else 
            isNewReadSegment = false;
      
        // if we've moved beyond current position
        if ( genomePosition > CurrentPosition ) {
            if ( op.Type == 'N' ) saveAlignment = false; // ignore alignment if REF_SKIP
            break;
        }
    }

    // save pileup position if flag is true
    if ( saveAlignment )
        CurrentPileupData.PileupAlignments.push_back( pileupAlignment );
}

// ---------------------------------------------
// PileupEngine implementation

PileupEngine::PileupEngine(void)
    : d( new PileupEnginePrivate )
{ }

PileupEngine::~PileupEngine(void) {
    delete d;
    d = 0;
}

bool PileupEngine::AddAlignment(const BamAlignment& al) { return d->AddAlignment(al); }
void PileupEngine::AddVisitor(PileupVisitor* visitor) { d->Visitors.push_back(visitor); }
void PileupEngine::Flush(void) { d->Flush(); }
