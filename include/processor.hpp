#pragma once
/**
 * FILE: processor.hpp
 * Created: Nov 06, 2017 Mon
 *
 * processing functions designed to run in threads
 */

// forward
class Readfeed;
struct Runopts;
struct Index;
struct Readstats;
class KeyValueDatabase;

void align(Readfeed& readfeed, Readstats& readstats, Index& index, KeyValueDatabase& kvdb, Runopts& opts);
