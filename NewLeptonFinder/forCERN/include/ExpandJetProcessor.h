#ifndef ExpandJetProcessor_h
#define ExpandJetProcessor_h 1

#include <string>

#include <marlin/Processor.h>
#include <lcio.h>

#include <EVENT/ReconstructedParticle.h>

using namespace lcio ;
using namespace marlin ;

class ExpandJetProcessor : public Processor {

	public:

		virtual Processor*  newProcessor() { return new ExpandJetProcessor ; }

		ExpandJetProcessor() ;

		virtual void init() ;
		virtual void processRunHeader( LCRunHeader* run ) ;
		virtual void processEvent( LCEvent * evt ) ; 
		virtual void check( LCEvent * evt ) ; 
		virtual void end() ;

	protected:

		/** Input collection */
		std::string _inputCollection;

		/** Output collections */
		std::string _outputCollection;

} ;

#endif

