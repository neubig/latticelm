/*
* Copyright 2010, Graham Neubig
* 
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
* 
*     http://www.apache.org/licenses/LICENSE-2.0
* 
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#ifndef LEXFST_H
#define LEXFST_H

#include "util.h"
#include <stdexcept>
#include <fst/vector-fst.h>

using namespace fst;

namespace latticelm {

template< class WordId, class CharId >
class LexFst : public VectorFst<StdArc> {

public:

    typedef typename StdArc::StateId StateId;

private:

    CharId numChars_;
    vector< vector<CharId> > words_;
    vector< string > symbols_;
    StateId homeState_, unkState_;
    string separator_;

public:

    LexFst() : VectorFst<StdArc>(), numChars_(0), words_(1), symbols_(), 
                homeState_(-1), unkState_(-1) {
    }
    
    // load symbols from a file and initialize the lexfst
    void load(const char* isymFile) {
        
        // initialize
        symbols_.push_back("<eps>");
        symbols_.push_back("<phi>");
        symbols_.push_back("x<unk>");
        symbols_.push_back("x</unk>");
        numChars_ = 2;

        // load from file
        ifstream in(isymFile);
        if(!in) {
            ostringstream err;
            err << "Could not find symbol file " << isymFile << endl;
            throw runtime_error(err.str());
        }
        string buff;
        // We need to make sure that the symbol input is appropriate for
        // latticelm
        vector<string> sanity(4);
        sanity[0] = "<eps>"; sanity[1] = "0";
        sanity[2] = "<phi>"; sanity[3] = "1";
        for(unsigned i = 0; i < 4; i++) {
            in >> buff;
            if(sanity[i] != buff)
                THROW_ERROR("The first two symbols in the symbol file must be \"<eps> 0\" and \"<phi> 1\"");
        }
        while(in >> buff) {
            // cerr << "Adding symbol " << buff << " as " << numChars_ << endl;
            numChars_++;
            symbols_.push_back("x"+buff);
            in >> buff;
        }
        symbols_.push_back("w<s>");
        initializeArcs();
    }

    void initializeArcs() {
        // initialize the states
        homeState_ = AddState(); 
        unkState_ = AddState(); 
        SetStart(homeState_);
        SetFinal(homeState_,0);
        AddArc(unkState_,StdArc(0,3,0,homeState_)); // end of unknown word
        for(unsigned i = 4; i < numChars_+2; i++) {
            StateId sid = AddState();
            AddArc(homeState_,StdArc(i-2,0,0,sid));
            AddArc(sid,StdArc(0,i,0,unkState_));
            AddArc(unkState_,StdArc(i-2,i,0,unkState_));
        }
    }

    vector<WordId> parseSample(const Fst<StdArc> & sample) {
        vector<WordId> ret;
        vector<CharId> charBuf;
        Fst<StdArc>::StateId sid = sample.Start();
        // continue until there are no more left
        while(true) {
            ArcIterator< Fst<StdArc> > ai(sample,sid);
            if(ai.Done())
                break;
            StdArc arc = ai.Value();
            // add known words
            if(arc.olabel >= numChars_+2) {
                if(charBuf.size() > 0)
                    throw std::runtime_error("Word with non-empty buffer (/unk required)");
                WordId wid = arc.olabel - numChars_-2;
                ret.push_back(wid);
            }
            // handle the end of unknown word symbol
            else if(arc.olabel == 3) {
                if(charBuf.size() == 0)
                    throw std::runtime_error("End of word symbol with empty buffer");
                charBuf.push_back(1);
                ret.push_back(addWord(charBuf));
                charBuf.resize(0);
            } 
            // handle unknown word characters
            else if(arc.olabel > 3)
                charBuf.push_back(arc.olabel-2);
            sid = arc.nextstate;
        }
        // clean up the remaining buffer
        if(charBuf.size() > 0)
            ret.push_back(addWord(charBuf));
        return ret;
    }

    
    WordId addWord(const vector<CharId> & word) {
        if(word.size() == 0)
            return 0;
        StateId sid = Start();
        CharId cid = 0;
        for(unsigned i = 0; i < word.size(); i++) {
            cid = word[i];
            if(cid == 1)
                break;
            ArcIterator< Fst<StdArc> > aiter(*this,sid);
            for(; !aiter.Done() && aiter.Value().ilabel != cid; aiter.Next());
            if(!aiter.Done())
                sid = aiter.Value().nextstate;
            else {
                StateId nextState = AddState();
                AddArc(sid, StdArc(cid,0,0,nextState));
                sid = nextState;
            }
        }
        // get the last value
        ArcIterator< Fst<StdArc> >  aiter(*this,sid);
        // find an arc that points to homeState_ && != cid+2
        for(; !aiter.Done() && 
            (aiter.Value().nextstate != homeState_ || aiter.Value().olabel == cid+2); 
            aiter.Next());
        // if the word is found, return its value
        if(!aiter.Done())
            return aiter.Value().olabel - numChars_ - 2;
        // otherwise, add the word to the dictionary
        WordId newId = words_.size()+numChars_+2;
        AddArc(sid, StdArc(0,newId,0,homeState_));
        words_.push_back(word);
        // add the symbol
        ostringstream oss;
        oss << 'w' << symbols_[word[0]+2].substr(1);
        for(unsigned i = 1; i < word.size()-1; i++) 
            oss << separator_ << symbols_[word[i]+2].substr(1);
        symbols_.push_back(oss.str());
        return words_.size()-1;
    }

    const vector< vector<CharId> > & getWords() { return words_; }
    const vector< string > & getSymbols() { return symbols_; }
    // get symbols that cannot be trimmed (character symbols + start/end symbols)
    vector< string > getPermSymbols() { 
        vector< string > ret(symbols_);
        ret.resize(numChars_+3);
        return ret;
    }
    void setPermSymbols(const vector< string > & symbols) { 
        symbols_ = symbols; 
        numChars_ = symbols_.size()-3;
    }
    void setSeparator(const string & separator) { separator_ = separator; }
    unsigned getNumChars() { return numChars_; }

};

}

#endif
