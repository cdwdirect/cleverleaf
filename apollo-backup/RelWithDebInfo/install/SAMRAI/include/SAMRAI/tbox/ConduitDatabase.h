/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   A database structure that stores (key,value) pairs in memory
 *
 ************************************************************************/

#ifndef included_tbox_ConduitDatabase
#define included_tbox_ConduitDatabase

#include "SAMRAI/SAMRAI_config.h"

#ifdef HAVE_CONDUIT

#include "SAMRAI/tbox/Database.h"

#include <list>
#include "conduit.hpp"

namespace SAMRAI {
namespace tbox {

/*!
 * Class ConduitDatabase stores (key,value) pairs in a hierarchical
 * database in memory using the Conduit Nodes for internal storage.
 * Each value may be another database, boolean, box, character,
 * complex, double, float, integer, or string.  Note that boxes are stored
 * using the DatabaseBox class that can store boxes of any dimension in the
 * same data structure.
 *
 * See the Database class documentation for a description of the
 * generic database interface.
 *
 * It is assumed that all processors will access the database in the same
 * manner.  Thus, all error messages are output to pout instead of perr.
 */

class ConduitDatabase:public Database
{
public:
   /*!
    * @brief Constructor creates an empty database with the specified name
    *
    * @param name  Name of the database
    */
   explicit ConduitDatabase(
      const std::string& name);

   /*!
    * @brief Constructor creates a database pointing to a Conduit Node
    *
    * @param name  Name of the database
    * @param node  Pointer to a Conduit Node--must be of Conduit's object
    *              data type.
    */
    ConduitDatabase(
       const std::string& name,
       conduit::Node* node);

   /*!
    * Destructor deallocates the data in the database.
    */
   virtual ~ConduitDatabase();

   /*!
    * @brief Fulfills abstract interface--not used in this implementation.
    *
    * An error will occur if this is called
    */
   virtual bool
   create(
      const std::string& name);

   /*!
    * @brief Fulfills abstract interface--not used in this implementation.
    *
    * An error will occur if this is called
    */
   virtual bool
   open(
      const std::string& name,
      const bool read_write_mode = false);

   /*!
    * @brief Fulfills abstract interface--not used in this implementation.
    *
    * An error will occur if this is called
    */
   virtual bool
   close();

   /**
    * Return string name of memory database object.
    */
   virtual std::string
   getName() const;

   /**
    * Return true if the specified key exists in the database and false
    * otherwise.
    */
   virtual bool
   keyExists(
      const std::string& key);

   /**
    * Return all keys in the database.
    */
   virtual std::vector<std::string>
   getAllKeys();

   /**
    * @brief Return the type of data associated with the key.
    *
    * If the key does not exist, then INVALID is returned
    *
    * @param key Key name in database.
    */
   virtual enum DataType
   getArrayType(
      const std::string& key);

   /**
    * Return the size of the array associated with the key.  If the key
    * does not exist, then zero is returned.
    */
   virtual size_t
   getArraySize(
      const std::string& key);

   /**
    * Return whether the specified key represents a database entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool
   isDatabase(
      const std::string& key);

   /**
    * Create a new database with the specified key name.  If the key already
    * exists in the database, then the old key record is deleted and the new
    * one is silently created in its place.
    */
   virtual std::shared_ptr<Database>
   putDatabase(
      const std::string& key);

   /**
    * Get the database with the specified key name.  If the specified
    * key does not exist in the database or it is not a database, then
    * an error message is printed and the program exits.
    */
   virtual std::shared_ptr<Database>
   getDatabase(
      const std::string& key);

   /**
    * Return whether the specified key represents a boolean entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool
   isBool(
      const std::string& key);

   /**
    * Create a boolean scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putBool(
      const std::string& key,
      const bool& data);

   /**
    * Create a boolean array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putBoolArray(
      const std::string& key,
      const bool * const data,
      const size_t nelements);

   /**
    * Get a boolean entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * boolean scalar, then an error message is printed and the program
    * exits.
    */
   virtual bool
   getBool(
      const std::string& key);

   /**
    * Get a boolean entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a boolean scalar,
    * then an error message is printed and the program exits.
    */
   virtual bool
   getBoolWithDefault(
      const std::string& key,
      const bool& defaultvalue);

   /**
    * Get a boolean entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a boolean vector, then an error message is printed and
    * the program exits.
    */
   virtual std::vector<bool>
   getBoolVector(
      const std::string& key);

   /**
    * Get a boolean entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a boolean array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void
   getBoolArray(
      const std::string& key,
      bool* data,
      const size_t nelements);

   /**
    * Return whether the specified key represents a box entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool
   isDatabaseBox(
      const std::string& key);

   /**
    * Create a box scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putDatabaseBox(
      const std::string& key,
      const DatabaseBox& data);

   /**
    * Create a box vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putDatabaseBoxVector(
      const std::string& key,
      const std::vector<DatabaseBox>& data);

   /**
    * Create a box array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putDatabaseBoxArray(
      const std::string& key,
      const DatabaseBox * const data,
      const size_t nelements);

   /**
    * Get a box entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * box scalar, then an error message is printed and the program
    * exits.
    */
   virtual DatabaseBox
   getDatabaseBox(
      const std::string& key);

   /**
    * Get a box entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a box scalar,
    * then an error message is printed and the program exits.
    */
   virtual DatabaseBox
   getDatabaseBoxWithDefault(
      const std::string& key,
      const DatabaseBox& defaultvalue);

   /**
    * Get a box entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a box vector, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual std::vector<DatabaseBox>
   getDatabaseBoxVector(
      const std::string& key);

   /**
    * Get a box entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a box array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void
   getDatabaseBoxArray(
      const std::string& key,
      DatabaseBox* data,
      const size_t nelements);

   /**
    * Return whether the specified key represents a character entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool
   isChar(
      const std::string& key);

   /**
    * Create a character scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putChar(
      const std::string& key,
      const char& data);

   /**
    * Create a character vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putCharVector(
      const std::string& key,
      const std::vector<char>& data);

   /**
    * Create a character array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putCharArray(
      const std::string& key,
      const char * const data,
      const size_t nelements);

   /**
    * Get a character entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * character scalar, then an error message is printed and the program
    * exits.
    */
   virtual char
   getChar(
      const std::string& key);

   /**
    * Get a character entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a character scalar,
    * then an error message is printed and the program exits.
    */
   virtual char
   getCharWithDefault(
      const std::string& key,
      const char& defaultvalue);

   /**
    * Get a character entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a character vector, then an error message is printed and
    * the program exits.
    */
   virtual std::vector<char>
   getCharVector(
      const std::string& key);

   /**
    * Get a character entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a character array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void
   getCharArray(
      const std::string& key,
      char* data,
      const size_t nelements);

   /**
    * Return whether the specified key represents a complex entry.  If
    * the key does not exist, then false is returned.  Complex values
    * may be promoted from integers, floats, or doubles.
    */
   virtual bool
   isComplex(
      const std::string& key);

   /**
    * Create a complex scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putComplex(
      const std::string& key,
      const dcomplex& data);

   /**
    * Create a complex vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putComplexVector(
      const std::string& key,
      const std::vector<dcomplex>& data);

   /**
    * Create a complex array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putComplexArray(
      const std::string& key,
      const dcomplex * const data,
      const size_t nelements);

   /**
    * Get a complex entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * complex scalar, then an error message is printed and the program
    * exits.  Complex values may be promoted from integers, floats, or
    * doubles.
    */
   virtual dcomplex
   getComplex(
      const std::string& key);

   /**
    * Get a complex entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a complex scalar,
    * then an error message is printed and the program exits.  Complex
    * values may be promoted from integers, floats, or doubles.
    */
   virtual dcomplex
   getComplexWithDefault(
      const std::string& key,
      const dcomplex& defaultvalue);

   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a complex vector, then an error message is printed and
    * the program exits.  Complex values may be promoted from integers,
    * floats, or doubles.
    */
   virtual std::vector<dcomplex>
   getComplexVector(
      const std::string& key);

   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a complex array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    * Complex values may be promoted from integers, floats, or doubles.
    */
   virtual void
   getComplexArray(
      const std::string& key,
      dcomplex* data,
      const size_t nelements);

   /**
    * Return whether the specified key represents a double entry.  If
    * the key does not exist, then false is returned.  Double values
    * may be promoted from integers or floats.
    */
   virtual bool
   isDouble(
      const std::string& key);

   /**
    * Create a double scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putDouble(
      const std::string& key,
      const double& data);

   /**
    * Create a double vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putDoubleVector(
      const std::string& key,
      const std::vector<double>& data);

   /**
    * Create a double array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putDoubleArray(
      const std::string& key,
      const double * const data,
      const size_t nelements);

   /**
    * Get a double entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * double scalar, then an error message is printed and the program
    * exits.  Double values may be promoted from integers or floats.
    */
   virtual double
   getDouble(
      const std::string& key);

   /**
    * Get a double entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a double scalar, then
    * an error message is printed and the program exits.  Double values may
    * be promoted from integers or floats.
    */
   virtual double
   getDoubleWithDefault(
      const std::string& key,
      const double& defaultvalue);

   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double vector, then an error message is printed and
    * the program exits.  Double values may be promoted from integers
    * or floats.
    */
   virtual std::vector<double>
   getDoubleVector(
      const std::string& key);

   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    * Double values may be promoted from integers or floats.
    */
   virtual void
   getDoubleArray(
      const std::string& key,
      double* data,
      const size_t nelements);

   /**
    * Return whether the specified key represents a float entry.  If
    * the key does not exist, then false is returned.  Float values
    * may be promoted from integers or silently truncated from doubles.
    */
   virtual bool
   isFloat(
      const std::string& key);

   /**
    * Create a float scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putFloat(
      const std::string& key,
      const float& data);

   /**
    * Create a float vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putFloatVector(
      const std::string& key,
      const std::vector<float>& data);

   /**
    * Create a float array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putFloatArray(
      const std::string& key,
      const float * const data,
      const size_t nelements);

   /**
    * Get a float entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * float scalar, then an error message is printed and the program
    * exits.  Float values may be promoted from integers or silently
    * truncated from doubles.
    */
   virtual float
   getFloat(
      const std::string& key);

   /**
    * Get a float entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a float scalar, then
    * an error message is printed and the program exits.  Float values may
    * be promoted from integers or silently truncated from doubles.
    */
   virtual float
   getFloatWithDefault(
      const std::string& key,
      const float& defaultvalue);

   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float vector, then an error message is printed and
    * the program exits.  Float values may be promoted from integers
    * or silently truncated from doubles.
    */
   virtual std::vector<float>
   getFloatVector(
      const std::string& key);

   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    * Float values may be promoted from integers or silently truncated
    * from doubles.
    */
   virtual void
   getFloatArray(
      const std::string& key,
      float* data,
      const size_t nelements);

   /**
    * Return whether the specified key represents an integer entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool
   isInteger(
      const std::string& key);

   /**
    * Create an integer scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putInteger(
      const std::string& key,
      const int& data);

   /**
    * Create an integer vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putIntegerVector(
      const std::string& key,
      const std::vector<int>& data);

   /**
    * Create an integer array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putIntegerArray(
      const std::string& key,
      const int * const data,
      const size_t nelements);

   /**
    * Get an integer entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * integer scalar, then an error message is printed and the program
    * exits.
    */
   virtual int
   getInteger(
      const std::string& key);

   /**
    * Get an integer entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not an integer scalar,
    * then an error message is printed and the program exits.
    */
   virtual int
   getIntegerWithDefault(
      const std::string& key,
      const int& defaultvalue);

   /**
    * Get an integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not an integer vector, then an error message is printed and
    * the program exits.
    */
   virtual std::vector<int>
   getIntegerVector(
      const std::string& key);

   /**
    * Get an integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not an integer array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void
   getIntegerArray(
      const std::string& key,
      int* data,
      const size_t nelements);

   /**
    * Return whether the specified key represents a std::string entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool
   isString(
      const std::string& key);

   /**
    * Create a string scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putString(
      const std::string& key,
      const std::string& data);

   /**
    * Create a string vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putStringVector(
      const std::string& key,
      const std::vector<std::string>& data);

   /**
    * Create a string array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putStringArray(
      const std::string& key,
      const std::string * const data,
      const size_t nelements);

   /**
    * Get a string entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * string scalar, then an error message is printed and the program
    * exits.
    */
   virtual std::string
   getString(
      const std::string& key);

   /**
    * Get a string entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a string scalar,
    * then an error message is printed and the program exits.
    */
   virtual std::string
   getStringWithDefault(
      const std::string& key,
      const std::string& defaultvalue);

   /**
    * Get a string entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a string vector, then an error message is printed and
    * the program exits.
    */
   virtual std::vector<std::string>
   getStringVector(
      const std::string& key);

   /**
    * Get a string entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a string array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void
   getStringArray(
      const std::string& key,
      std::string* data,
      const size_t nelements);

   /**
    * @brief Returns the name of this database.
    *
    * The name for the root of the database is the name supplied when creating
    * it.  Names for nested databases are the keyname of the database.
    */
   virtual std::string
   getName();

   /**
    * Print the current database to the specified output stream.  After
    * each key, print whether that key came from the a file and was
    * used, came from the file but was not used (unused),
    * or came from a default key value (default).  If no output stream
    * is specified, then data is written to stream pout.
    *
    * NOTE:  under the g++ compiler libraries, printClassData has a
    * maximum output of 4096 characters per line.
    */
   virtual void
   printClassData(
      std::ostream& os = pout);

   /**
    * Print the database keys that were not used to the specified output
    * stream.
    */
   void
   printUnusedKeys(
      std::ostream& os = pout) const
   {
      printDatabase(os, 0, PRINT_UNUSED);
   }

   /**
    * Print the database keys that were set via default calls to the specified
    * output stream.
    */
   void
   printDefaultKeys(
      std::ostream& os = pout) const
   {
      printDatabase(os, 0, PRINT_DEFAULT);
   }

private:

   ConduitDatabase();                            // not implemented
   ConduitDatabase(
      const ConduitDatabase&);                   // not implemented
   ConduitDatabase&
   operator = (
      const ConduitDatabase&);                   // not implemented

   /*
    * Private utility routines for managing the database
    */
   bool
   deleteKeyIfFound(
      const std::string& key);

   conduit::Node&
   getChildNodeOrExit(
      const std::string& key);


   void
   setConduitDataTypes();

   static void
   indentStream(
      std::ostream& os,
      const long indent)
   {
      for (int i = 0; i < indent; ++i) {
         os << " ";
      }
   }

   void
   printDatabase(
      std::ostream& os,
      const int indent,
      const int toprint) const;

   /*
    * Private data members - name and a list of (key,value) pairs
    */
   std::string d_database_name;

   conduit::Node* d_node;

   std::map<std::string, std::shared_ptr<ConduitDatabase> > d_child_dbs;
   std::map<std::string, enum Database::DataType> d_types;

   static const int PRINT_DEFAULT;
   static const int PRINT_INPUT;
   static const int PRINT_UNUSED;
   static const int SSTREAM_BUFFER;
};

}
}

#endif

#endif
