#ifndef KQUANTITY_H
#define KQUANTITY_H

//STL headers
#include <string>
#include <sstream>

//custom headers
#include "Analysis/KCode/KMap.h"

using namespace std;

//------------------------------------------
//generic base class for a quantity
class KQuantity {
	public:
		//constructor
		KQuantity() : name(""), localOpt(new OptionMap()) { }
		KQuantity(string name_, OptionMap* localOpt_) : name(name_), localOpt(localOpt_ ? localOpt_ : new OptionMap()) { }
		//destructor
		virtual ~KQuantity() {}
		
		//accessors
		virtual string GetName() { return name; }
		virtual string LegName() { return legname; }
		virtual string PrintName() { return printname; }
		virtual string AxisName() { return axisname; }
		virtual string LegVsName() { return legvsname; }
		virtual string PrintVsName() { return printvsname; }
		virtual string CutName() { return cutname; }
		virtual string GetDesc(bool print, bool varied){
			if(print){
				if(varied) return printvsname;
				else return printname;
			}
			else {
				if(varied) return legvsname;
				else return legname;
			}
		}
		virtual double GetVal() { return 0.; }
		virtual double GetErrLow() { return 0.; }
		virtual double GetErrHigh() { return 0.; }
		OptionMap* GetLocalOpt() { return localOpt; }
		void SetLocalOpt(OptionMap* opt) { localOpt = opt ? opt : new OptionMap(); }
		
		//operators (for sorting, etc.)
		bool operator<(const KQuantity& rhs){
			return this->name < rhs.name;
		}
		bool operator==(const KQuantity& rhs){
			return Compare(rhs);
		}
		bool operator!=(const KQuantity& rhs){
			return !(*this == rhs);
		}
		virtual bool Compare(const KQuantity& rhs){
			return this->name == rhs.name;
		}
		
	protected:
		//member variables
		string name;
		OptionMap* localOpt;
		string legname, printname, axisname, legvsname, printvsname, cutname;
};

#endif