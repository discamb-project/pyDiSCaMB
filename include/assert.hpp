// Custom assert to throw Python AssertionError instead of abort

#include <exception>
#include <sstream>

class AssertionError : public std::exception {
    public:
        AssertionError(const char *expr, const char *file, int line){
            std::stringstream ss;
            ss << "AssertionError: " << expr << "\nThrown from line " << line << " in file " << file;
            mMessage = ss.str();
        };
        const char * what() const noexcept {
            return mMessage.c_str();
        };

    private:
        std::string mMessage;
};

#ifdef assert
    #undef assert
#endif
#define assert(expr) (static_cast <bool> (expr) ? void (0) : throw AssertionError(#expr, __FILE__, __LINE__))
