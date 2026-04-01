#ifndef TOKENIZER_HH
#define TOKENIZER_HH

#include <string>

template <class ContainerT>
void tokenize(const std::string_view &str, ContainerT &tokens,
              const std::string_view &delimiters = " \n\t", bool trimEmpty = true) {
    std::string::size_type lastPos = 0, length = str.length();

    using value_type = typename ContainerT::value_type;
    using size_type = typename ContainerT::size_type;

    while(lastPos < length + 1) {
        std::string::size_type pos = str.find_first_of(delimiters, lastPos);
        if(pos == std::string::npos) pos = length;

        if(pos != lastPos || !trimEmpty) {
            tokens.push_back(
                value_type(str.data() + lastPos, static_cast<size_type>(pos) - lastPos));
        }
        lastPos = pos + 1;
    }
}

#endif
