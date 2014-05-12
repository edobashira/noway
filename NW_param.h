// This may look like C code, but it is really -*- C++ -*-
#include <assert.h>
#include <memory.h>
#include <stdio.h>
#include <string.h>
#include "NW_type.h"


// Parameter tables consist of a set of parameter entries
struct PARAM {
  PARAM(const char* _name) : name(_name) { next = NULL; data = NULL; size = 0; }

  PARAM* next;  // pointer to next parameter entry in list
  const char* name;         // name of parameter
  int size;           // number of strings under this parameter name
  CHAR** data;        // list of strings
  void resize(int new_size);  // change size of list of strings
};


enum PARAM_MODE { OVERRIDE, UNDERRIDE, APPEND };

class UNCHECKED_PARAM_TABLE {
 public:
  UNCHECKED_PARAM_TABLE(const char* file_name = NULL); 

  virtual ~UNCHECKED_PARAM_TABLE(); 

  // Parse a C-style command line argument list
  void parse_cmd_line(INT argc, const char** argv, PARAM_MODE mode = OVERRIDE);

  // Parse parameters from a file
  virtual void parse_file(const char* filename, PARAM_MODE mode = OVERRIDE);

  // Parse a string of parameter(s)
  BOOL parse_line(CHAR* string, PARAM_MODE mode = OVERRIDE);

  // Get the value of a parameter (or return default_value if not found)
  // Each parameter can have a list of strings; index is used to
  // select a string from this list.
  BOOL get_bool(const char* name);

  // return number of tokens in parameter (= max(index) - 1)
  INT get_size(const char* name);

  INT get_int(const char* name, INT default_value = 0, INT index = 0);


  float get_float(const char* name, float default_value = 0.0, INT index = 0);

  const char* get_string(const char* name, const char* default_value = "", INT index = 0);
  CHAR* copy_string(const char* name, const char* default_value = "", INT index = 0);

  INT get_enum(const char* name, INT default_value, 
               const char* *enum_name_table, INT index = 0);

  // get array of ints or floats
  INT* get_int_array(const char* name);
  float* get_float_array(const char* name);
  INT* get_range_array(const char* name, INT size, INT included_in_range = 1, 
		       INT excluded_from_range = 0);

  INT* get_mapping(const char* name, INT *map_size, INT not_mapped = -1);

//  void set_path_prefix(NAMED*);  // find() returns longest matching path


  // print and log functions
  void print(FILE* file = stdout, STRING first_char = "");  //DCA 12/MAR/96
  void print_entry(FILE* file, PARAM* entry, STRING first_char = ""); //DCA 12/MAR/96

  void set_default_log_file(FILE* log_file) { default_log = log_file; }

  // Set default parameter name (used before the first -name arg)
  void set_default_param_name(const char* name) 
    { default_param_name = name; }

 protected:
  virtual PARAM *put(const char* name);   // put a parameter in the table
  virtual PARAM *find(const char* name);  // find existing parameter(else NULL)
  void del(const char* name);        // remove a name from table

  // Set parameter name [ index ] to value.
  void put(const char* name, int index, const char* value);

  // Start up a parse through a list of parameters.
  // (if no initial -parameter_name, then assign values to parameter)
  void start_parse() { current_param = strdup(default_param_name); }

  // Parse a single word from the input file or command argument list.
  void parse_token(const char* token, PARAM_MODE mode);

 private:
  PARAM *head, *tail;       // head and tail of PARAM list

  // Parameter name to which current values are being added.
  char* current_param;

  // Initial parameter name
  const char* default_param_name;
  FILE* default_log;

 protected:
  CHAR* file_name;
};

class PARAM_TABLE : public UNCHECKED_PARAM_TABLE {
 public:
  PARAM_TABLE(const char* file_name = NULL, const char** param_table = NULL);

  PARAM* put(const char* name);

  PARAM* find(const char* name);

  void help_message();

 private:
  INT check_param;
  UNCHECKED_PARAM_TABLE master_parameter_list;

  INT get_num_parameters(const char** p); 
};


extern const char* itoa(int);  // convert int to string


