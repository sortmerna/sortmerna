#pragma once
/**
 * FILE: processor.hpp
 * Created: Nov 06, 2017 Mon
 *
 * processing functions designed to run in threads
 */
#include <string>
#include <vector>

// forward
class Read;
class Readfeed;
struct Runopts;
struct Index;
class References;
class Output;
struct Readstats;
class Refstats;
class KeyValueDatabase;

void align(Readfeed& readfeed, Readstats& readstats, Index& index, KeyValueDatabase& kvdb, Runopts& opts);
void align2(int id, Readfeed& readfeed, Readstats& readstats, Index& index, References& refs, Refstats& refstats, KeyValueDatabase& kvdb, Runopts& opts);
void writeSummary(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts); /* post-alignment tasks like calculating statistics */
void writeSummary2(int id, Readfeed& readfeed, Runopts& opts, References& refs, Readstats& readstats, Refstats& refstats, KeyValueDatabase& kvdb);
void writeSummary3(Read& read, Readstats& readstats, Refstats& refstats, References& refs, Runopts& opts);
void writeReports(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts);
void report(int id, Readfeed& readfeed, Runopts& opts, References& refs, Refstats& refstats, Output& output, KeyValueDatabase& kvdb);
