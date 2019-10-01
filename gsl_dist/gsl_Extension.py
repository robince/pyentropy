#
# author: Achim Gaedke
# created: May 2001
# file: pygsl/gsl_dist/gsl_extension.py
# $Id: gsl_Extension.py,v 1.4 2008/02/09 16:42:58 schnizer Exp $
#
# module for gsl extensions compilation

from distutils.core import setup, Extension
from distutils.errors import DistutilsModuleError, DistutilsExecError
import os
import os.path
import re
import string
import types
import imp
from sys import argv,version_info



from gsl_dist.array_includes import array_include_dirs


# steel --gsl-prefix from option list
gsl_prefix_option=None
gsl_prefix_option_pattern=re.compile("--gsl-prefix=(.+)")
pos=0
while pos<len(argv):
	gsl_prefix_match=gsl_prefix_option_pattern.match(argv[pos])
	if gsl_prefix_match:
		gsl_prefix_option=gsl_prefix_match.group(1)
		gsl_prefix_option.strip()
		argv[pos:pos+1]=[]
		break
	pos+=1

# Extension class adapted for gsl
class _gsl_Location:
	"""
	Wrapper for the location of the gsl library.

	On unix one can run gsl-config to find the locations. On other systems
	one has to revert to other ways to find the configuration.
	"""
	def __init__(self):
		self.prefix = None
		self.cflags = None
		self.libs   = None
		self.version = None
		self.swig = None
		
	def get_gsl_prefix(self):
		assert(self.prefix != None)
		return self.prefix
			
	def get_gsl_cflags(self):
		assert(self.cflags != None)
		return self.cflags
	
	def get_gsl_libs(self):
		#print self.libs
		assert(self.libs != None)
		return self.libs
	
	
	def get_gsl_version(self):
		assert(self.version != None)
		return self.version
	
	def _split_version(self, version):
		if version[-1] == '+':
			version = version[:-1]
		return re.split('\.',version)

	def get_swig(self):
		assert(self.swig)
		return self.swig
	
class _gsl_Location_gsl_config(_gsl_Location):
	"""
	Call gsl_config to find the location of gsl
	"""
	def __init__(self):
		_gsl_Location.__init__(self)
		gsl_prefix = None
		if gsl_prefix is not None:
			self.gsl_config_tool=os.path.join(gsl_prefix,"bin","gsl-config")
		elif gsl_prefix_option is not None:
			self.gsl_config_tool=os.path.join(gsl_prefix_option,"bin","gsl-config")
		else:
			self.gsl_config_tool="gsl-config"
			
		self.prefix = self.get_gsl_info('--prefix').strip()
		self.cflags = self.get_gsl_info('--cflags').strip()
		self.libs   = self.get_gsl_info('--libs').strip()
		self.version = self._split_version(self.get_gsl_info('--version').strip())[:2]
		
		# I am running on swig. I suppose that swig is in the path
		self.swig = "swig"
		try:
			self.swig = os.environ["SWIG"]
		except KeyError:
			pass

		
	def get_gsl_info(self, arguments):
		"""
		executes gsl-config with given arguments
		"""
		gsl_command=os.popen(self.gsl_config_tool+' '+arguments)
		gsl_output=gsl_command.readline()
		gsl_command.close()
		if not gsl_output:
			raise DistutilsExecError("could not start %s"%self.gsl_config_tool)
		return  gsl_output

	
class _gsl_Location_file(_gsl_Location):
	def __init__(self):
		_gsl_Location.__init__(self)
		try:
			import gsl_site
		except ImportError as des:
			msg = "I do not know how to run gsl-config \n"+\
				  "on this system. Therefore you must provide me with the information\n" +\
				  "where to find the GSL library. I could not import `gsl_site'.\n" +\
				  "Reason: %s. Copy gsl_site_example.py to gsl_site.py.\n"+\
				  "Edit the variables in that file to reflect your installation."
			raise DistutilsExecError(msg % des)
		
		self.prefix  = gsl_site.prefix 
		self.cflags  = gsl_site.cflags 
		self.libs	= gsl_site.libs   
		self.swig	= gsl_site.swig
		self.version = self._split_version(gsl_site.version)

if os.name == 'posix':
	gsl_Location = _gsl_Location_gsl_config()
else:
	gsl_Location = _gsl_Location_file()


class gsl_Extension(Extension):
	"""
	for gsl needs
	"""
	def __init__(self, name, sources,
		include_dirs=None,
		define_macros=None,
		undef_macros=None,
		library_dirs=None,
		libraries=None,
		runtime_library_dirs=None,
		extra_objects=None,
		extra_compile_args=None,
		extra_link_args=None,
		export_symbols=None,
		gsl_prefix=None,
		gsl_min_version=None,
		python_min_version=None):
	
		# get real prefix
		self.gsl_prefix=self.get_gsl_prefix()
		
		gsl_major_version, gsl_minor_version = self.get_gsl_version()
		# check gsl version
		if gsl_min_version is not None and \
			not self.check_gsl_version(gsl_min_version):
				raise DistutilsExecError(\
				  "min gsl version %s required"%repr(gsl_min_version))

		# check python version
		if python_min_version is not None and \
			not self.check_python_version(python_min_version):
				raise DistutilsExecError( \
				  "min python version %s required"%repr(python_min_version))

		# prepend include directory
		if include_dirs is None: include_dirs=[]
		include_dirs.append('Include')
		include_dirs.append('.')
		include_dirs[0:0]=[os.path.join(self.gsl_prefix,'include')]
		include_dirs= include_dirs + array_include_dirs

		# prepend library directory
		if library_dirs is None: library_dirs=[]
		library_dirs[0:0] = [os.path.join(self.gsl_prefix,'lib')]

		# prepare lib list
		# split optionlist and strip blanks from each option
		gsl_lib_list=map(string.strip,self.get_gsl_libs().split())

		# filter options with -l
		not_lib_opt=lambda a:a[:2]=="-l"
		gsl_lib_list=filter(not_lib_opt,gsl_lib_list)
		# cut -l
		only_lib_name=lambda a:a[2:]
		gsl_lib_list=map(only_lib_name,gsl_lib_list)

		if libraries is None: libraries=[]
		#libraries.append('pygsl')
		libraries.extend(gsl_lib_list)

		# test if Numeric module is available
		if define_macros is None:
			define_macros=[]
		try:
			imp.find_module("Numeric")
			define_macros = define_macros + [("NUMERIC",1),]
		except ImportError:
			define_macros = define_macros + [("NUMERIC",0), ]
		if undef_macros == None:
			undef_macros = []
		if 'NDEBUG' not in undef_macros:
			undef_macros.append('NDEBUG')
		tmp = map(lambda x: x[0], define_macros)
		if "PYGSL_GSL_MAJOR_VERSION" not in tmp:
			define_macros = define_macros + [("PYGSL_GSL_MAJOR_VERSION", gsl_major_version),]
		if "PYGSL_GSL_MINOR_VERSION" not in tmp:
			#define_macros.append(("PYGSL_GSL_MINOR_VERSION", gsl_minor_version))
			define_macros = define_macros + [("PYGSL_GSL_MINOR_VERSION", gsl_minor_version),]
		
		Extension.__init__(self, name, sources,
			include_dirs,
			define_macros,
			undef_macros,
			library_dirs,
			libraries,
			runtime_library_dirs,
			extra_objects,
			extra_compile_args,
			extra_link_args,
			export_symbols
			)

	def check_version(self, required_version, this_version):
		min_length=min(len(required_version),len(this_version))
		for pos in range(min_length):
			this_type=type(required_version[pos])
			if  this_type== types.StringType:
				if required_version[pos]>this_version[pos]: return 0
			elif this_type == types.IntType:
				if required_version[pos]>int(this_version[pos]): return 0
			else:
				raise DistutilsExecError("incorrect version specification")
				# problematic: 0.9 < 0.9.0, but assures that 0.9.1 > 0.9
				if len(required_version)>len(this_version): return 0
		return 1


	def check_gsl_version(self, version_array):
		return self.check_version(version_array,self.get_gsl_version())

	def check_python_version(self, version_array):
		return self.check_version(version_array,version_info)
			
	# get gsl-prefix option
	def get_gsl_info(self, arguments):
		"""
		executes gsl-config with given arguments
		"""
		gsl_command=os.popen(self.gsl_config_tool+' '+arguments)
		gsl_output=gsl_command.readline()
		gsl_command.close()
		if not gsl_output:
			raise DistutilsExecError("could not start %s"%self.gsl_config_tool)
		return  gsl_output

	def get_gsl_prefix(self,):
		return gsl_Location.get_gsl_prefix()

	def get_gsl_cflags(self):
		return gsl_Location.get_gsl_cflags()

	def get_gsl_libs(self):
		return gsl_Location.get_gsl_libs()

	def get_gsl_version(self):
		return gsl_Location.get_gsl_version()

