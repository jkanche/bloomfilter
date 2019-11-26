package main 

import (
	"fmt"
	"math/rand"
	"math"
	"os"
	"bufio"
	"strings"
	"strconv"
	"encoding/json"
	"time"
	"encoding/gob"
	"github.com/spaolacci/murmur3"
)

const CACHE_SIZE = 512

// bitvec structure
type bitvec struct {
	Bits []uint8
}

// create a bit vector
func makeBitvec(input []bool) *bitvec {
	size := (len(input) / 8)

	if len(input) % 8 > 0 {
		size += 1
	}

	// fmt.Println("size ", size)
	bvec := new(bitvec)
	bvec.Bits = make([]uint8, size)

	val := 0

	if input[0] == true {
		val += int(math.Pow(float64(2), float64(7)))
	}
	// fmt.Println("val ", val)
	// fmt.Println("len(input) ", len(input))

	for i := 1; i < len(input); i++ {

		// fmt.Println("i ", i)

		if i/8 == 0 {
			// fmt.Println("i/8 == 0 ", i/8)

			bvec.Bits[i/8] = uint8(val)
			val = 0
		}

		mod := 7 - (i % 8)

		// fmt.Println("mod ", mod)

		if input[i] == true {
			val += int(math.Pow(float64(2), float64(mod)))
		}
	}

	return bvec
}

// use bitwise AND to get bit at a position
func (bv bitvec) getBitAt(i int) uint8 {
	block := i / 8
	mod := i % 8

	bitcheck := uint8(math.Pow(float64(2), float64(mod)))

	if bv.Bits[block] & bitcheck == 0 {
		return 0
	}
	
	return 1
}

func (bv bitvec) setBitAt(i int) {
	block := i / 8
	mod := i % 8

	bitcheck := uint8(math.Pow(float64(2), float64(mod)))

	if bv.Bits[block] & bitcheck == 0 {
		bv.Bits[block] += bitcheck
	}
}

type bloomfilter struct {
	Bitvec *bitvec
	Bitsize int
	fpr float64
	Hash_size int
	Hash_seeds []uint32
}

func makeBloomFilter(bitsize int, fpr float64, k int) *bloomfilter {
	bloom := new(bloomfilter)
	size := (bitsize / 8)

	if bitsize % 8 > 0 {
		size += 1
	}

	bloom.Bitvec = new(bitvec)
	bloom.Bitvec.Bits = make([]uint8, size)
	bloom.Hash_size = k
	bloom.Hash_seeds = make([]uint32, k)
	bloom.Bitsize = bitsize

	for i := 0; i < k; i++ {
		bloom.Hash_seeds[i] = uint32(rand.Intn(100000))
	}

	// fmt.Println("bloom", bloom)
	return bloom
}

func (bloom bloomfilter) insert(key string) {

	// choose block
	block := int((int(murmur3.Sum32WithSeed([]byte(key), bloom.Hash_seeds[0])) % bloom.Bitsize) / CACHE_SIZE)

	for j := 1; j < bloom.Hash_size; j++ {
		// fmt.Println("j, hash_seed", j, bloom.Hash_seeds[j])
		bitval := int((int(murmur3.Sum32WithSeed([]byte(key), bloom.Hash_seeds[j])) % bloom.Bitsize) % CACHE_SIZE)
		// fmt.Println("bit val", bitval)
		bitpos := block * CACHE_SIZE + bitval
		// fmt.Println(bitpos)
		bloom.Bitvec.setBitAt(bitpos)
	}
}

func (bloom bloomfilter) query(key string) string {

	var sum uint8 = 0

	// choose block
	block := int((int(murmur3.Sum32WithSeed([]byte(key), bloom.Hash_seeds[0])) % bloom.Bitsize) / CACHE_SIZE)

	for j := 1; j < bloom.Hash_size; j++ {
		// fmt.Println("j, hash_seed", j, bloom.Hash_seeds[j])
		bitval := int((int(murmur3.Sum32WithSeed([]byte(key), bloom.Hash_seeds[j])) % bloom.Bitsize) % CACHE_SIZE)
		// fmt.Println("bit val", bitval)
		bitpos := block * CACHE_SIZE + bitval
		sum += bloom.Bitvec.getBitAt(bitpos)
	}

	if int(sum) == bloom.Hash_size - 1 {
		return "Y"
	}

	return "N"
}

func parse_opts(args []string) map[string]string {

	if len(args) % 2 != 0 {
		fmt.Println("incorrect params!")
		panic("incorrect params")
	}

	argmap := make(map[string]string)
	
	for i := 0; i < len(args); i++ {
		argmap[args[i]] = args[i+1]
		i++
	}

	return argmap
}

func (bloom bloomfilter) save_bf(fileName string) {
	f, err := os.Create(fileName)
	if err != nil {
		panic(err)
	}

	enc := gob.NewEncoder(f)
	err = enc.Encode(bloom)
	if err != nil {
		panic(err)
	}
	f.Close()
}

// loads bloom from disk
func load_bf(fileName string) *bloomfilter {
	f, err := os.Open(fileName)
	if err != nil {
		panic(err)
	}
	dec := gob.NewDecoder(f)

	var bloom bloomfilter
	dec.Decode(&bloom)
	f.Close()

	return &bloom
}

func run_tests() {
	
	var presenttimes, abscenttimes, mixedtimes, absfpr, mixfpr [4][5]float64

	// i -> size of the data, starts at 10K 
	// j-> fpr starts at 1
	for i := 1; i < 5; i++ {
		for j := 1; j < 6; j++ {
			nKeys := 1000 * int(math.Pow(float64(10), float64(i)))
			fpr := float64(j) * float64(0.05)
			modk := int(nKeys/1000)

			bitvec_size := int(float64(-1.0) * float64(nKeys) * math.Log(fpr) / math.Pow(math.Log(2), 2))

			size := (bitvec_size / CACHE_SIZE)
			if bitvec_size % CACHE_SIZE > 0 {
				size += 1
			}

			bitvec_size = size * CACHE_SIZE
			hash_func_k := int(float64(-1.0) * math.Log2(fpr))
			
			file, err := os.Open("../../kmer.shuf")
			if err != nil {
				panic(err)
			}

			fmt.Println("nKeys, fpr, bitvec_size, hash_func_k", nKeys, fpr, bitvec_size, hash_func_k)

			bloom := makeBloomFilter(bitvec_size, fpr, hash_func_k)

			scanner := bufio.NewScanner(file)
			count, abscount := 0, 0
			var present, abscent [1000]string

			for scanner.Scan() {
				key := strings.TrimSpace(scanner.Text())

				if abscount == 1000 {
					break
				}

				if count < nKeys {
					bloom.insert(key)

					// add every 100th word to known list
					if count % modk == 0 {
						// fmt.Println(count)
						present[int(count / modk)] = key
					}

					count++
				} else {
					abscent[abscount] = key
					abscount++
				}
				// count++
				// if count > nKeys {
				// 	abscent
				// }
			}

			// fmt.Println("present", len(present))
			// fmt.Println(present)
			// fmt.Println("abscent", len(abscent))
			// fmt.Println(abscent)

			start := time.Now()
			pfpr := 0
			for k := 0; k < len(present); k++ {
				res := bloom.query(present[k])
				if res == "N" {
					pfpr++
				}
			}
			end := time.Now().Sub(start)
			presenttimes[i-1][j-1] = float64(end.Nanoseconds())
			fmt.Println("pfpr, ", pfpr)

			start = time.Now()
			afpr := 0
			for k := 0; k < len(abscent); k++ {
				res := bloom.query(abscent[k])
				if res == "Y" {
					afpr++
				}
			}
			end = time.Now().Sub(start)
			abscenttimes[i-1][j-1] = float64(end.Nanoseconds())
			absfpr[i-1][j-1] = float64(afpr)

			start = time.Now()
			mfpr := 0
			for k := 0; k < 500; k++ {
				bloom.query(present[k])
			}
			for k := 0; k < 500; k++ {
				res := bloom.query(abscent[k])
				if res == "Y" {
					mfpr++
				}
			}
			end = time.Now().Sub(start)
			mixedtimes[i-1][j-1] = float64(end.Nanoseconds())
			mixfpr[i-1][j-1] = float64(mfpr)
		}
	}

	// presenttimes, abscenttimes, mixedtimes, absfpr, mixfpr
	f, _ := os.Create("presenttimes.json")
	enc := json.NewEncoder(f)
	_ = enc.Encode(presenttimes)
	f.Close()

	f, _ = os.Create("abscenttimes.json")
	enc = json.NewEncoder(f)
	_ = enc.Encode(abscenttimes)
	f.Close()

	f, _ = os.Create("mixedtimes.json")
	enc = json.NewEncoder(f)
	_ = enc.Encode(mixedtimes)
	f.Close()

	f, _ = os.Create("absfpr.json")
	enc = json.NewEncoder(f)
	_ = enc.Encode(absfpr)
	f.Close()

	f, _ = os.Create("mixfpr.json")
	enc = json.NewEncoder(f)
	_ = enc.Encode(mixfpr)
	f.Close()
}


func main() {

	// first argument is the name of the binary
	args := os.Args[1:]

	// parse the operation to perform
	switch args[0] {
		case "build":
			argmap := parse_opts(args[1:])

			fmt.Println("argmap", argmap)
			
			fpr, err := strconv.ParseFloat(argmap["-f"], 64)
			nKeys, err := strconv.ParseInt(argmap["-n"], 10, 0)

			fmt.Println("fpr, nkeys", fpr, nKeys)
			// calculate size 
			bitvec_size := int(float64(-1.0) * float64(nKeys) * math.Log(fpr) / math.Pow(math.Log(2), 2))
			hash_func_k := int(float64(-1.0) * math.Log2(fpr))

			fmt.Println("bitvec_size, hash_func_k", bitvec_size, hash_func_k)

			bloom := makeBloomFilter(bitvec_size, fpr, hash_func_k)
			
			fmt.Println("bloom created", bloom)

			file, err := os.Open(argmap["-k"])

			if err != nil {
				panic(err)
			}

			scanner := bufio.NewScanner(file)
			for scanner.Scan() {
				bloom.insert(scanner.Text())
			}

			if err := scanner.Err(); err != nil {
				panic(err)
			}

			bloom.save_bf(argmap["-o"])

		case "query":
			argmap := parse_opts(args[1:])

			bloom := load_bf(argmap["-i"])

			file, err := os.Open(argmap["-q"])
			if err != nil {
				panic(err)
			}

			scanner := bufio.NewScanner(file)
			for scanner.Scan() {
				fmt.Println(bloom.query(scanner.Text()))
			}

			if err := scanner.Err(); err != nil {
				panic(err)
			}
		case "runtests":
			run_tests()
		default:
			fmt.Println("unrecognized command!")
	}
}