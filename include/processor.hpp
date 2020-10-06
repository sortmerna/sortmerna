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

void align(Runopts& opts, Readstats& readstats, Output& output, Index& index, KeyValueDatabase& kvdb);
void align2(int id, Readfeed& readfeed, Index& index, References& refs, Readstats& readstats, Refstats& refstats, KeyValueDatabase& kvdb, Runopts& opts);
void postProcess(Runopts& opts, Readstats& readstats, Output& output, KeyValueDatabase& kvdb); /* post-alignment tasks like calculating statistics */
void postProcess2(int id, Readfeed& readfeed, Runopts& opts, References& refs, Readstats& readstats, Refstats& refstats, KeyValueDatabase& kvdb);
void computeStats(Read& read, Readstats& readstats, Refstats& refstats, References& refs, Runopts& opts);
void generateReports(Runopts& opts, Readstats& readstats, Output& output, KeyValueDatabase& kvdb);
void reportsJob(Readfeed& readfeed, Runopts& opts, References& refs, Refstats& refstats, Output& output, KeyValueDatabase& kvdb);
