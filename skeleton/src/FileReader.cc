#include <FileReader.hh>
#include "Debug.hh"

void FileReader::registerIntParameter(const std::string &key, int init)
{
  IntParameter[key] = init;
}

void FileReader::registerRealParameter(const std::string &key, real init)
{
  RealParameter[key] = init;
}

void FileReader::registerStringParameter(const std::string &key, const std::string &init)
{
  StringParameter[key] = init;
}

void FileReader::setParameter(const std::string &key, const std::string &in)
{
 	 std::string message = {"StringParameter not registered: "};
 	 message += key; 	 
 	 CHECK_MSG(StringParameter.count(key) != 0,message);
   StringParameter[key] = in;
}

void FileReader::setParameter(const std::string &key, real in)
{
 	 std::string message = {"RealParameter not registered: "};
 	 message += key;
 	 CHECK_MSG(RealParameter.count(key) != 0,message); 	 
   RealParameter[key] = in;
}

void FileReader::setParameter(const std::string &key, int in)
{
 	 std::string message = {"IntParameter not registered: "};
 	 message += key;
 	 CHECK_MSG(IntParameter.count(key) != 0,message); 	 
   IntParameter[key] = in;
}


bool FileReader::readFile(const std::string &name)
{
        std::ifstream ifs(name.c_str());
		std::string buffer;
		std::string key;
		std::string value;
    long long ii = 0;
		if(ifs.is_open()){
			while (ifs.good())
			{
                ++ii;
				getline(ifs,buffer);
				if(FileReader::has_only_spaces(buffer))
				{
					continue;
				}
				size_t strlength = buffer.find("#") - buffer.find("\n");
				if (strlength!=0)
				{
					buffer.erase(buffer.find("#"));
				}
				std::istringstream iss (buffer);
				iss >> std::skipws >> key  >> std::skipws >> value;
                std::string err_msg = "Malformed input file! Could not read parameter: ";
                err_msg += key;
                err_msg += "\n In line: ";
                err_msg += std::to_string(ii);
                if(key.length()==0) continue;
                CHECK_MSG(value.length()!=0,err_msg);

				if(IntParameter.count(key)!=0)
				{
                    std::string err_msg1 = "Parameter is of wrong type (registered as int): ";
                    err_msg1 += key;
                    err_msg1 += "\n In line: ";
                    err_msg1 += std::to_string(ii);
                    CHECK_MSG(string_only_digits(value),err_msg1);
					int valuei = stoi(value);
					FileReader::setParameter(key,valuei);
				}
				else if(RealParameter.count(key)!=0)
				{
                    std::string err_msg1 = "Parameter is of wrong type (registered as real): ";
                    err_msg1 += key;
                    err_msg1 += "\n In line: ";
                    err_msg1 += std::to_string(ii);
                    CHECK_MSG(string_only_reals(value),err_msg1);
					real valuei = stof(value);
					FileReader::setParameter(key,valuei);
				}
				else if(StringParameter.count(key)!=0)
				{
					FileReader::setParameter(key,value);
				}
				key.clear();
				value.clear();	  
			}
			ifs.close();
			return true;
		}
		else
		return false;
}

bool FileReader::has_only_spaces(const std::string& str) {
	if(str.find_first_not_of (' ') == str.npos||str.find_first_not_of ('\t') == str.npos)
		return true;
	else
		return false;
}

bool FileReader::string_only_digits(const std::string& str){
    return str.find_first_not_of("0123456789-") == std::string::npos;
}

bool FileReader::string_only_reals(const std::string& str){
    return str.find_first_not_of("0123456789.e-") == std::string::npos;
}

void FileReader::printParameters() const
{
    std::cout << "Printing out all integer parameters: \n";
	for (std::map<std::string, int>::const_iterator i = IntParameter.begin(); i != IntParameter.end(); ++i)
	{
		std::cout << i->first << " " << i->second << std::endl;
	}
	std::cout << std::endl;
    std::cout << "Printing out all real parameters: \n";
	for (std::map<std::string, real>::const_iterator i = RealParameter.begin(); i != RealParameter.end(); ++i)
	{
		std::cout << i->first << " " << i->second << std::endl;
	}
	std::cout << std::endl;
    std::cout << "Printing out all string parameters: \n";
	for (std::map<std::string, std::string>::const_iterator i = StringParameter.begin(); i != StringParameter.end(); ++i)
	{
		std::cout << i->first << " " << i->second << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Finished printing all parameters \n\n";
}


