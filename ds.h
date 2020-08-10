#ifndef DS_h
#define DS_h

class Site {
public:
	int chr, pos;
	Site() {}
	Site(int Chr, int Pos) {
		chr = Chr;
		pos = Pos;
	}
	bool operator > (const Site &rhs ) const {
		return ( chr > rhs.chr || (chr==rhs.chr && pos>rhs.pos) );
	}
	bool operator < (const Site &rhs ) const {
		return ( chr < rhs.chr || (chr==rhs.chr && pos<rhs.pos) );
	}
	bool operator == (const Site &rhs ) const {
		return ( chr==rhs.chr && pos==rhs.pos );
	}
};

class Interval {
public:
	int chr, start, end;
	Interval() {}
	Interval(int Chr, int Start, int End) {
		chr = Chr;
		start = Start;
		end = End;
	}
	bool startsInside(const Interval &window) const {
		return (chr==window.chr && window.start <= start && start <= window.end);
	}
	bool spansInto(const Interval &window) const {
		return (chr==window.chr && start < window.start && end >= window.start);
	}
	bool endsInside(const Interval &window) const {
		return (chr==window.chr && window.start <= end && end <= window.end);
	}
	bool overlapsWith(const Interval &window) const {
		return (this->startsInside(window) || this->spansInto(window));
	}
	bool operator > (const Interval &interval) const {
		return (chr > interval.chr || (chr==interval.chr && start>interval.start));
	}
	bool operator < (const Interval &interval) const {
		return (chr < interval.chr || (chr==interval.chr && start<interval.start));
	}
	bool operator == (const Interval &interval) const {
		return (chr==interval.chr && start==interval.chr && end==interval.end);
	}
};

#endif

