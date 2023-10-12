package main

import (
	"bufio"
	"errors"
	"flag"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"sort"
	"strconv"
	"strings"
	"time"

	"gonum.org/v1/gonum/stat"
)

var n int
var file_name string
var N int
var M int
var TS int
var P int
var r float64
var two_opt_rate float64
var T float64

var arg_n *int
var arg_file_name *string
var arg_N *int
var arg_M *int
var arg_r *float64
var arg_T *float64
var arg_two_opt_rate *float64
var arg_best_tour_length *float64
var arg_TS *int
var arg_init_alg *string
var arg_report_period *int
var arg_action *string

type InitPopulationFunc func(int) (error, int, []Solution)
type DoCrossover func() error

type Point struct {
	x float64
	y float64
}

type Pair struct {
	x int
	y int
}

type Solution struct {
	h          float64
	chromosome []int
	ord        []int
}

type Edge struct {
	to     int
	length float64
}

func (e Edge) String() string {
	return fmt.Sprintf("%d: %f", e.to, e.length)
}

type Node struct {
	id    int
	alpha int
	pos   int
}

type ByLength []Edge

func (a ByLength) Len() int           { return len(a) }
func (a ByLength) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a ByLength) Less(i, j int) bool { return a[i].length < a[j].length }

var population []Solution
var best_solution int
var best_h = 1e10
var edge_weight [][]float64
var sorted_node [][]int
var count_evol_from_crossover = 0
var count_evol_from_mutation = 0
var count_evol_from_two_opt = 0
var report_period = 1000
var best_tour []int
var best_tour_lenght float64
var opt_tour_avai bool
var best_tour_mask []int
var population_edge_mask [][]int
var combined_population []Solution
var is_nearest_rank bool

func read_node_file(file_name string) []Point {
	node_list := make([]Point, n)
	nodefile, err := os.Open(file_name + ".tsp")
	if err != nil {
		log.Fatal(err)
	}
	defer nodefile.Close()
	scanner := bufio.NewScanner(nodefile)

	start := false
	for scanner.Scan() {
		line := scanner.Text()
		tokens := strings.Fields(line)
		if start {
			if tokens[0] == "EOF" {
				break
			} else if len(tokens) == 3 {
				i, _ := strconv.Atoi(tokens[0])
				i -= 1
				if i >= n {
					break
				}
				x, _ := strconv.ParseFloat(tokens[1], 64)
				y, _ := strconv.ParseFloat(tokens[2], 64)
				node_list[i] = Point{x, y}
				// fmt.Println(node_list[i])
			}
		} else if tokens[0] == "NODE_COORD_SECTION" || tokens[0] == "DISPLAY_DATA_SECTION" {
			start = true
		} else if tokens[0] == "EDGE_WEIGHT_SECTION" {
			return nil
		}
	}

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}
	return node_list
}

func read_edge_weight(file_name string) [][]float64 {
	edge_weight := make([][]float64, n)
	for i, _ := range edge_weight {
		edge_weight[i] = make([]float64, n)
	}

	nodefile, err := os.Open(file_name + ".tsp")
	if err != nil {
		log.Fatal(err)
	}
	defer nodefile.Close()

	scanner := bufio.NewScanner(nodefile)

	start := false
	matrixType := 0
	i := 0
	j := 0
	for scanner.Scan() {
		line := scanner.Text()
		tokens := strings.Fields(line)
		if start {
			if tokens[0] == "EOF" || tokens[0] == "DISPLAY_DATA_SECTION" {
				break
			} else {
				if matrixType == 0 {
					for k := 0; k < len(tokens); k++ {
						w, e := strconv.Atoi(tokens[k])
						if e == nil {
							edge_weight[i][j] = float64(w)
							j++
							if j >= n {
								j = 0
								i++
							}
						}
					}
				}
				if matrixType == 1 {
					for k := 0; k < len(tokens); k++ {
						w, e := strconv.Atoi(tokens[k])
						if e == nil {
							j++
							if j >= n {
								i++
								j = i + 1
							}
							edge_weight[i][j] = float64(w)
							edge_weight[j][i] = float64(w)
						}
					}
				}
				if matrixType == 2 {
					for k := 0; k < len(tokens); k++ {
						w, e := strconv.Atoi(tokens[k])
						if e == nil {
							edge_weight[i][j] = float64(w)
							edge_weight[j][i] = float64(w)
							j++
							if j >= n {
								i++
								j = i
							}
						}
					}
				}
				if matrixType == 3 {
					// fmt.Println("tokens", tokens)
					for k := 0; k < len(tokens); k++ {
						w, e := strconv.Atoi(tokens[k])
						if e == nil {
							if i < n && j < n {
								edge_weight[i][j] = float64(w)
								edge_weight[j][i] = float64(w)
								// } else {
								// 	fmt.Println("tokens", tokens)
							}
							// fmt.Println(i, j, tokens[k])
							// for i := 0; i < len(edge_weight); i++ {
							// 	fmt.Println(edge_weight[i])
							// }
							j++
							if j > i {
								j = 0
								i++
							}
						}
					}
				}
			}
		} else if tokens[0] == "EDGE_WEIGHT_SECTION" {
			start = true
			fmt.Println("EDGE_WEIGHT_SECTION", tokens)
		} else if tokens[0] == "EDGE_WEIGHT_FORMAT" || tokens[0] == "EDGE_WEIGHT_FORMAT:" {
			if tokens[len(tokens)-1] == "UPPER_ROW" {
				matrixType = 1
			} else if tokens[len(tokens)-1] == "UPPER_DIAG_ROW" {
				matrixType = 2
			} else if tokens[len(tokens)-1] == "LOWER_DIAG_ROW" {
				matrixType = 3
			} else {
				matrixType = 0
			}
		} else {
			// fmt.Println("tokens", tokens)
		}
	}
	// fmt.Println("matrixType", matrixType)
	// for i := 0; i < len(edge_weight); i++ {
	// 	fmt.Println(edge_weight[i])
	// }

	return edge_weight
}

func read_tour_file(file_name string) []int {
	tour := make([]int, n)
	tourfile, err := os.Open(file_name + ".opt.tour")
	if err != nil {
		fmt.Println(err)
		return nil
	}
	defer tourfile.Close()
	scanner := bufio.NewScanner(tourfile)

	start := false
	i := -1
	for scanner.Scan() {
		line := scanner.Text()
		tokens := strings.Fields(line)
		// fmt.Println("tokens",tokens)
		if start {
			if tokens[0] == "-1" {
				break
			} else {
				for k := 0; k < len(tokens); k++ {
					j, _ := strconv.Atoi(tokens[k])
					if i >= n {
						break
					}
					tour[i] = j - 1
					i += 1
				}
				// fmt.Println(nodeList[i])
			}
		} else if tokens[0] == "TOUR_SECTION" {
			start = true
			i = 0
		}
	}

	return tour
}

func write_history(file_name string, history []float64) {
	f, _ := os.Create(file_name)
	w := bufio.NewWriter(f)
	for i := 0; i < len(history); i++ {
		w.WriteString(strconv.Itoa(int(history[i])) + "\n")
	}
	w.Flush()
	f.Close()
}

func write_tour(file_name string, tour []int) {
	f, _ := os.Create(file_name)
	w := bufio.NewWriter(f)
	for i := 0; i < len(tour); i++ {
		w.WriteString(strconv.Itoa(int(tour[i])) + "\n")
	}
	w.Flush()
	f.Close()
}

func writePopulationH(file_name string) {
	f, _ := os.Create(file_name)
	w := bufio.NewWriter(f)
	for i := 0; i < len(population); i++ {
		w.WriteString(strconv.Itoa(int(population[i].h)) + "\n")
	}
	w.Flush()
	f.Close()
}

func calc_nearest() {
	if sorted_node == nil {
		sorted_node = make([][]int, n)
	}
	for i := 0; i < n; i++ {
		sorted_node[i] = make([]int, n-1)

		edges := make(ByLength, n)
		for j := 0; j < n; j++ {
			edges[j] = Edge{j, edge_weight[i][j]}
		}
		// fmt.Println(edges[i])
		sort.Sort(edges)
		ind := 0
		for j := 0; j < n; j++ {
			if edges[j].to != i {
				sorted_node[i][ind] = edges[j].to
				ind++
			}
		}
	}
	is_nearest_rank = false
}

func calc_nearest_rank() {
	node_rank := make([][]int, n)
	for i := 0; i < n; i++ {
		node_rank[i] = make([]int, n)
		edges := make(ByLength, n)
		for j := 0; j < n; j++ {
			edges[j] = Edge{j, edge_weight[i][j]}
		}
		// fmt.Println(edges[i])
		sort.Sort(edges)
		ind := 0
		for j := 0; j < n; j++ {
			if edges[j].to != i {
				node_rank[i][edges[j].to] = ind
				ind++
			}
		}
	}
	if sorted_node == nil {
		sorted_node = make([][]int, n)
	}
	for i := 0; i < n; i++ {
		sorted_node[i] = make([]int, n-1)
		edges := make(ByLength, n)
		for j := 0; j < n; j++ {
			rank := float64(node_rank[i][j] + node_rank[j][i])
			edges[j] = Edge{j, rank}
		}
		sort.Sort(ByLength(edges))
		ind := 0
		for j := 0; j < n; j++ {
			if edges[j].to != i {
				sorted_node[i][ind] = edges[j].to
				ind++
			}
		}
	}
	is_nearest_rank = true
}

func calc_edge_weight(node_list []Point) [][]float64 {
	edge_weight := make([][]float64, n)
	for i, _ := range edge_weight {
		edge_weight[i] = make([]float64, n)
	}
	// fmt.Println(edge_weight[0][0])
	// fmt.Println(edge_weight[1][1])
	// fmt.Println(edge_weight[2][2])
	// fmt.Println(edge_weight[3][3])
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			// fmt.Println(i, j)
			pi := node_list[i]
			pj := node_list[j]
			edge_weight[i][j] = math.Sqrt(math.Pow(pi.x-pj.x, 2) + math.Pow(pi.y-pj.y, 2))
			edge_weight[j][i] = math.Sqrt(math.Pow(pi.x-pj.x, 2) + math.Pow(pi.y-pj.y, 2))
		}
	}
	return edge_weight
}

func calc_tour_lenght(tour []int) (tour_lenght float64) {
	// tour_lenght := 0
	for i := 0; i < n-1; i++ {
		tour_lenght += edge_weight[tour[i]][tour[i+1]]
	}
	tour_lenght += edge_weight[tour[n-1]][tour[0]]
	return tour_lenght
}

func calc_tour_lenght2(tour []int) (tour_lenght float64) {
	// tour_lenght := 0
	for i := 0; i < n-2; i++ {
		tour_lenght += edge_weight[tour[i]][tour[i+1]]
	}
	tour_lenght += edge_weight[0][tour[0]]
	tour_lenght += edge_weight[tour[n-2]][0]
	return tour_lenght
}

func get_ord(chromosome []int) []int {
	result := make([]int, n)
	for i := 0; i < n; i++ {
		result[chromosome[i]] = i
	}
	return result
}

func check(chromosome []int) bool {
	mask := make([]bool, n)
	for i := 0; i < n; i++ {
		mask[chromosome[i]] = true
	}
	for i := 0; i < n; i++ {
		if !mask[i] {
			fmt.Println("lack", i)
			return false
		}
	}
	return true
}

func gsx2(solparent1 Solution, solparent2 Solution, p int) ([]int, []int) {
	taken := make([]bool, n)
	child := make([]int, n)
	ord1 := solparent1.ord
	ord2 := solparent2.ord
	parent1 := solparent1.chromosome
	parent2 := solparent2.chromosome
	idx1 := (ord1[p] + 1) % n
	idx2 := (ord2[p] - 1 + n) % n
	idx := 1
	idx3 := 1
	idx4 := n - 1
	done := false
	taken[p] = true
	child[0] = p
	i := 0
	// fmt.Println("====================================================")
	for idx < n {
		i++

		// fmt.Println(child)
		if i > n*3 {
			fmt.Println("Danger !!!!")
			fmt.Println(parent1)
			fmt.Println(parent2)
			fmt.Println(p)
			fmt.Println(child)
			break
		}
		if done {
			if !taken[parent2[idx2]] {
				child[idx4] = parent2[idx2]
				taken[parent2[idx2]] = true
				idx += 1
				idx4 = (idx4 - 1 + n) % n
			}
			idx2 = (idx2 - 1 + n) % n
		} else if parent1[idx1] == parent2[idx2] {
			done = true
		} else if taken[parent1[idx1]] && (!taken[parent2[idx2]]) {
			done = true
			child[idx4] = parent2[idx2]
			taken[parent2[idx2]] = true
			idx1 = (idx1 + 1) % n
			idx2 = (idx2 - 1 + n) % n
			idx4 = (idx4 - 1 + n) % n
			idx += 1
		} else if (!taken[parent1[idx1]]) && taken[parent2[idx2]] {
			done = true
			child[idx3] = parent1[idx1]
			taken[parent1[idx1]] = true
			idx1 = (idx1 + 1) % n
			idx2 = (idx2 - 1 + n) % n
			idx3 = (idx3 + 1) % n
			idx += 1
		} else if taken[parent1[idx1]] && taken[parent2[idx2]] {
			done = true
		} else {
			child[idx3] = parent1[idx1]
			child[idx4] = parent2[idx2]
			taken[parent1[idx1]] = true
			taken[parent2[idx2]] = true
			idx += 2
			idx1 = (idx1 + 1) % n
			idx2 = (idx2 - 1 + n) % n
			idx3 = (idx3 + 1) % n
			idx4 = (idx4 - 1 + n) % n
		}
	}
	// if !check(child) {
	// 	fmt.Println("Toang !!!!")
	// 	fmt.Println(parent1)
	// 	fmt.Println(parent2)
	// 	fmt.Println(p)
	// 	fmt.Println(child)
	// }
	return child, get_ord(child)
}

func gsx3(solparent1 Solution, solparent2 Solution, p int) ([]int, []int) {
	taken := make([]bool, n)
	child := make([]int, n)
	ord1 := solparent1.ord
	ord2 := solparent2.ord
	parent1 := solparent1.chromosome
	parent2 := solparent2.chromosome
	child_done := 0
	candidate := p
	P = 1
	for i := 0; i < 2*n; i++ {
		taken[candidate] = true
		child[child_done] = candidate
		child_done++
		last_visit := candidate
		if i < P {
			candidate = parent2[(ord2[last_visit]+1)%n]
			if taken[candidate] {
				candidate = parent1[(ord1[last_visit]+1)%n]
			}
			if taken[candidate] {
				candidate = -1
			}
		} else {
			candidate = parent1[(ord1[last_visit]+1)%n]
			if taken[candidate] {
				candidate = parent2[(ord2[last_visit]+1)%n]
			}
			if taken[candidate] {
				candidate = -1
			}
		}
		if candidate == -1 {
			for j := 0; j < n-1; j++ {
				candidate = sorted_node[last_visit][j]
				if !taken[candidate] {
					break
				}
			}
		}
		if child_done >= n {
			break
		}
	}

	// fmt.Println("p", p)
	// fmt.Println("parent1", parent1)
	// fmt.Println("parent2", parent2)
	// fmt.Println("child", child)
	if !check(child) {
		fmt.Println("Toang !!!!")
		// fmt.Println(parent1)
		// fmt.Println(parent2)
		// fmt.Println(p)
		// fmt.Println(child)
	}
	return child, get_ord(child)
}

func pmx1(parent1 []int, parent2 []int, p1 int, p2 int) []int {
	child1 := make([]int, n)
	taken := make([]int, n)
	mask := make([]int, n)
	for j := p1; j <= p2; j++ {
		mask[parent2[j]] = 1
		taken[parent2[j]] = 1
	}
	for j := 0; j < p1; j++ {
		taken[parent1[j]] = 1
	}
	for j := p2 + 1; j < n; j++ {
		taken[parent1[j]] = 1
	}
	for i := 0; i < n; i++ {
		if i < p1 || i > p2 {
			if mask[parent1[i]] == 1 {
				for j := p1; j <= p2; j++ {
					if taken[parent1[j]] == 0 {
						child1[i] = parent1[j]
						taken[child1[i]] = 1
						break
					}
				}
			} else {
				child1[i] = parent1[i]
				taken[child1[i]] = 1
			}
		} else {
			child1[i] = parent2[i]
			taken[child1[i]] = 1
		}
	}
	return child1
}

func pmx(solparent1 Solution, solparent2 Solution) ([]int, []int) {
	parent1 := solparent1.chromosome
	parent2 := solparent2.chromosome

	p1 := rand.Intn(n)
	p2 := rand.Intn(n - 1)
	if p2 >= p1 {
		p2++
	} else {
		p2, p1 = p1, p2
	}
	child1 := pmx1(parent1, parent2, p1, p2)
	child2 := pmx1(parent2, parent1, p1, p2)

	if !check(child1) || !check(child2) {
		fmt.Println("Toang !!!!")
		fmt.Println(parent1)
		fmt.Println(child1)
		fmt.Println(parent2)
		fmt.Println(child2)
		fmt.Println(parent1[p1 : p2+1])
		fmt.Println(parent2[p1 : p2+1])
	}
	return child1, child2
}

func ANNmutation(parentf Solution) ([]int, []int) {
	parent := parentf.chromosome
	ord := parentf.ord
	mask := make([]int, n)
	fchoice := rand.Intn(n - 1)
	mask[fchoice] = 1
	left := parent[(ord[fchoice]+1+n)%n]
	right := parent[(ord[fchoice]-1+n)%n]
	mask[left] = 1
	mask[right] = 1
	prop := make([]float64, n-1)
	sum := float64(0)
	for j := 0; j < n-1; j++ {
		candidate := sorted_node[fchoice][j]
		prop[j] = float64(1-mask[candidate]) * math.Exp(-float64(edge_weight[fchoice][candidate])/100)
		sum += prop[j]
	}
	for j := 0; j < n-1; j++ {
		prop[j] = prop[j] / sum
	}
	choice := rand_choice_w_prop(n-1, prop)
	schoice := sorted_node[fchoice][choice]

	child := make([]int, n)
	copy(child, parent)
	picked := ord[left]
	iter1 := picked
	child[iter1] = schoice
	iter1 = (iter1 + 1) % n
	iter2 := ord[schoice]
	iter2 = (iter2 - 1 + n) % n
	for parent[iter2] != fchoice {
		child[iter1] = parent[iter2]
		iter1 = (iter1 + 1) % n
		iter2 = (iter2 - 1 + n) % n
	}
	return child, get_ord(child)
}

func rand_choice_w_prop(n int, prop []float64) int {
	choice := -1
	rand_choice := rand.Float64()
	// fmt.Println("rand_choice", rand_choice)
	cumulation := prop[0]
	for i := 1; i < n; i++ {
		// fmt.Println(i, "cumulation", cumulation)
		if cumulation >= rand_choice {
			choice = i - 1
			return choice
		}
		cumulation += prop[i]
	}
	return n - 1
}

func random_initilization(n_solution int) (error, int, []Solution) {
	best_h = 1e10
	best_pos := -1
	population := make([]Solution, n_solution)
	for c := 0; c < n_solution; c++ {
		gene := rand.Perm(n)

		if !check(gene) {
			return errors.New("!check(gene)"), -1, nil
		}
		h := calc_tour_lenght(gene)
		population[c].chromosome = gene
		population[c].ord = get_ord(gene)
		population[c].h = h
		if h < best_h {
			best_pos = c
			best_h = h
		}
	}
	return nil, best_pos, population
}

func nn_initilization(n_solution int) (error, int, []Solution) {
	best_h = 1e10
	best_pos := -1
	population := make([]Solution, n_solution)
	for c := 0; c < n_solution; c++ {
		gene := make([]int, n)
		mask := make([]int, n)
		fchoice := rand.Intn(n)
		gene[0] = fchoice
		mask[fchoice] = 1
		for i := 1; i < n; i++ {
			available_nearest := make([]int, n)
			available_count := 0
			for j := 0; j < n-1; j++ {
				k := sorted_node[gene[i-1]][j]
				// if k == 133
				if mask[k] == 0 {
					if available_count == 0 || (available_count > 0 && (int(edge_weight[gene[i-1]][k]) == int(edge_weight[gene[i-1]][available_nearest[0]]))) {
						available_nearest[available_count] = k
						available_count += 1
					} else {
						break
					}
				}
			}
			if available_count == 0 {
				fmt.Println("detected")
			}
			kk := rand.Intn(available_count)
			gene[i] = available_nearest[kk]
			mask[gene[i]] = 1
		}

		if !check(gene) {
			return errors.New("!check(gene)"), -1, nil
		}
		h := calc_tour_lenght(gene)
		population[c].chromosome = gene
		population[c].ord = get_ord(gene)
		population[c].h = h
		if h < best_h {
			best_pos = c
			best_h = h
			// fmt.Println("h", h, " best_h ", best_h)
		}
	}
	return nil, best_pos, population
}

func multi_nn_initilization(n_solution int) (error, int, []Solution) {
	best_h = 1e10
	best_pos := -1
	population := make([]Solution, n_solution)
	for c := 0; c < n_solution; c++ {
		gene := make([]int, n)
		mask1 := make([]int, n)
		for i := 0; i < n; i++ {
			mask1[i] = 2
		}
		cluster := make([]int, n)
		for i := 0; i < n; i++ {
			cluster[i] = -1
		}
		adj := make([][]bool, n)
		for i := 0; i < n; i++ {
			adj[i] = make([]bool, n)
			for j := 0; j < n; j++ {
				adj[i][j] = true
			}
		}

		for i := 0; i < n-1; i++ {
			prop := make([]float64, n)
			for j := 0; j < n; j++ {
				prop[j] = float64(mask1[j])
			}
			sum_prop := float64(0)
			for j := 0; j < n; j++ {
				sum_prop += prop[j]
			}
			for j := 0; j < n; j++ {
				prop[j] /= sum_prop
			}
			fchoice := rand_choice_w_prop(n, prop)
			// fmt.Println("fchoice", fchoice)
			if fchoice == -1 {
				fmt.Println("Toang", n, prop)
			}
			if cluster[fchoice] == -1 {
				cluster[fchoice] = fchoice
			}
			mask1[fchoice] -= 1

			mask := make([]bool, n)
			for j := 0; j < n; j++ {
				mask[j] = (mask1[j] > 0) && adj[fchoice][j] && ((cluster[j] - cluster[fchoice]) != 0)
			}
			mask[fchoice] = false
			schoice := -1

			available_nearest := make([]int, n)
			available_count := 0
			for j := 0; j < n-1; j++ {
				k := sorted_node[fchoice][j]
				// if k == 133
				if mask[k] {
					if available_count == 0 || (available_count > 0 && (int(edge_weight[fchoice][k]) == int(edge_weight[fchoice][available_nearest[0]]))) {
						available_nearest[available_count] = k
						available_count += 1
					} else {
						break
					}
				}
			}
			if available_count == 0 {
				return errors.New("available_count == 0"), -1, nil
			}
			kk := rand.Intn(available_count)
			schoice = available_nearest[kk]

			mask1[schoice] -= 1
			if cluster[schoice] == -1 {
				cluster[schoice] = cluster[fchoice]
			} else {
				cls := cluster[schoice]
				for j := 0; j < n; j++ {
					if cluster[j] == cls {
						cluster[j] = cluster[fchoice]
					}
				}
			}
			adj[fchoice][schoice] = false
			adj[schoice][fchoice] = false
		}
		for i := 0; i < n-1; i++ {
			if mask1[i] > 0 {
				for j := i + 1; j < n; j++ {
					if mask1[j] > 0 {
						adj[i][j] = false
						adj[j][i] = false
						mask1[i] = 0
						mask1[j] = 0
						break
					}
				}
				break
			}
		}

		gene[0] = 0
		for j := 0; j < n; j++ {
			if !adj[0][j] {
				gene[1] = j
				break
			}
		}
		for i := 2; i < n; i++ {
			last := gene[i-1]
			for j := 0; j < n; j++ {
				if (!adj[last][j]) && j != gene[i-2] {
					gene[i] = j
					break
				}
			}
		}

		if !check(gene) {
			return errors.New("!check(gene)"), -1, nil
		}
		h := calc_tour_lenght(gene)
		population[c].chromosome = gene
		population[c].ord = get_ord(gene)
		population[c].h = h
		if h < best_h {
			best_h = h
			best_pos = c
		}
	}
	return nil, best_pos, population
}

func rnn_initilization(n_solution int) (error, int, []Solution) {
	best_h = 1e10
	best_pos := -1
	population := make([]Solution, n_solution)
	for c := 0; c < n_solution; c++ {
		gene := make([]int, n)
		mask := make([]int, n)
		fchoice := rand.Intn(n)
		// fchoice := c
		// fchoice := fchoices[c]
		// fmt.Println("==========================================")
		// fmt.Println("fchoice", fchoice)
		gene[0] = fchoice
		mask[fchoice] = 1
		for i := 1; i < n; i++ {
			// fmt.Println("==========================================")
			for j := 0; j < n; j++ {
				if j >= n-1 {
					fmt.Println("==========================================")
					fmt.Println("detected")
					fmt.Println("fchoice", fchoice)
					fmt.Println("i=", i, "j=", j)
					fmt.Println(gene)
					fmt.Println("gene[i-1]", gene[i-1])
					// fmt.Println("gene=", gene)
					sort.Ints(gene)
					fmt.Println("len(sorted_gene)=", len(gene))
					fmt.Println("mask=", mask)
					for k := 1; k < len(gene); k++ {
						if gene[k] != gene[k-1]+1 {
							fmt.Println("k=", k, "sorted_gene[k]", gene[k])
							fmt.Println(gene[k-1 : k+10])
							// break
						}
					}
				} else {
					// if k == 133
					k := sorted_node[gene[i-1]][j]
					// fmt.Println("k", k, mask[k])
					if mask[k] == 0 {
						gene[i] = k
						mask[gene[i]] = 1
						break
					}
				}
			}
		}

		if !check(gene) {
			return errors.New("!check(gene)"), -1, nil
		}
		// population[c] = encode(gene)
		// raw_population[c] = gene
		h := calc_tour_lenght(gene)
		population[c].chromosome = gene
		population[c].ord = get_ord(gene)
		population[c].h = h
		// hpopulation[c] = int(h)
		// fmt.Println("=================================================================")
		// fmt.Println("gene", c, gene)
		// fmt.Println("population[c]", population[c])
		// fmt.Println("h", h)
		if h < best_h {
			// fmt.Println("replace", h, best_solution.h)
			best_h = h
			best_pos = c
		}
	}
	return nil, best_pos, population
}

func multi_rnn_initilization(n_solution int) (error, int, []Solution) {
	best_h = 1e10
	best_pos := -1
	population := make([]Solution, n_solution)
	for c := 0; c < n_solution; c++ {
		// fmt.Println("=================================================================")
		gene := make([]int, n)
		mask1 := make([]int, n)
		for i := 0; i < n; i++ {
			mask1[i] = 2
		}
		cluster := make([]int, n)
		for i := 0; i < n; i++ {
			cluster[i] = -1
		}
		adj := make([][]bool, n)
		for i := 0; i < n; i++ {
			adj[i] = make([]bool, n)
			for j := 0; j < n; j++ {
				adj[i][j] = true
			}
		}

		for i := 0; i < n-1; i++ {
			prop := make([]float64, n)
			for j := 0; j < n; j++ {
				prop[j] = float64(mask1[j])
			}
			sum_prop := float64(0)
			for j := 0; j < n; j++ {
				sum_prop += prop[j]
			}
			for j := 0; j < n; j++ {
				prop[j] /= sum_prop
			}
			fchoice := rand_choice_w_prop(n, prop)
			// fmt.Println("fchoice", fchoice)
			if fchoice == -1 {
				fmt.Println("Toang", n, prop)
			}
			if cluster[fchoice] == -1 {
				cluster[fchoice] = fchoice
			}
			mask1[fchoice] -= 1

			mask := make([]bool, n)
			for j := 0; j < n; j++ {
				mask[j] = (mask1[j] > 0) && adj[fchoice][j] && ((cluster[j] - cluster[fchoice]) != 0)
			}
			mask[fchoice] = false
			schoice := -1
			for j := 0; j < n; j++ {
				k := sorted_node[fchoice][j]
				if mask[k] {
					schoice = k
					break
				}
			}
			mask1[schoice] -= 1
			if cluster[schoice] == -1 {
				cluster[schoice] = cluster[fchoice]
			} else {
				cls := cluster[schoice]
				for j := 0; j < n; j++ {
					if cluster[j] == cls {
						cluster[j] = cluster[fchoice]
					}
				}
			}
			adj[fchoice][schoice] = false
			adj[schoice][fchoice] = false
			// fmt.Println("fchoice", fchoice, "schoice", schoice)
			// fmt.Println("mask1", mask1)
			// fmt.Println("prop", prop)
		}
		for i := 0; i < n-1; i++ {
			if mask1[i] > 0 {
				for j := i + 1; j < n; j++ {
					if mask1[j] > 0 {
						adj[i][j] = false
						adj[j][i] = false
						mask1[i] = 0
						mask1[j] = 0
						break
					}
				}
				break
			}
		}

		gene[0] = 0
		for j := 0; j < n; j++ {
			if !adj[0][j] {
				gene[1] = j
				break
			}
		}
		for i := 2; i < n; i++ {
			last := gene[i-1]
			for j := 0; j < n; j++ {
				if (!adj[last][j]) && j != gene[i-2] {
					gene[i] = j
					break
				}
			}
		}

		h := calc_tour_lenght(gene)
		population[c].chromosome = gene
		population[c].ord = get_ord(gene)
		population[c].h = h
		// hpopulation[c] = int(h)
		// fmt.Println("=================================================================")
		// fmt.Println("gene", c, gene)
		// fmt.Println("population[c]", population[c])
		// fmt.Println("h", h)
		if h < best_h {
			// fmt.Println("replace", h, best_solution.h)
			best_h = h
			best_pos = c
		}
	}
	return nil, best_pos, population
}

func do_crossover() error {
	p1 := rand.Intn(N)
	// ================================================
	p2 := rand.Intn(N - 1)
	if p2 >= p1 {
		p2 += 1
	}

	p := rand.Intn(n)
	var child1, ord1, child2, ord2 []int

	child1, ord1 = gsx3(population[p1], population[p2], p)
	child2, ord2 = gsx3(population[p2], population[p1], p)

	hchild1 := calc_tour_lenght(child1)
	hchild2 := calc_tour_lenght(child2)

	hchild := hchild2
	if hchild1 < hchild2 {
		hchild = hchild1
	}
	if hchild < population[p1].h && hchild < population[p2].h {
		count_evol_from_crossover++
		copy(population[p1].chromosome, child1)
		copy(population[p1].ord, ord1)
		population[p1].h = hchild1
		if hchild1 < best_h {
			best_solution = p1
			best_h = hchild1
			best_tour = child1
		}

		copy(population[p2].chromosome, child2)
		copy(population[p2].ord, ord2)
		population[p2].h = hchild2
		if hchild2 < best_h {
			best_solution = p2
			best_h = hchild2
			best_tour = child2
		}
	}

	return nil
}

func do_crossover_pmx() error {
	p1 := rand.Intn(N)
	p2 := rand.Intn(N - 1)
	if p2 >= p1 {
		p2 += 1
	}
	parent1 := population[p1]
	parent2 := population[p2]

	child1, child2 := pmx(parent1, parent2)
	ord1 := get_ord(child1)
	ord2 := get_ord(child2)

	// fmt.Println(child2, rawchild2)
	hchild1 := calc_tour_lenght(child1)
	hchild2 := calc_tour_lenght(child2)

	hchild := hchild2
	if hchild1 < hchild2 {
		hchild = hchild1
	}
	if hchild < population[p1].h && hchild < population[p2].h {
		count_evol_from_crossover++
		copy(population[p1].chromosome, child1)
		copy(population[p1].ord, ord1)
		population[p1].h = hchild1
		if hchild1 < best_h {
			best_solution = p1
			best_h = hchild1
		}

		copy(population[p2].chromosome, child2)
		copy(population[p2].ord, ord2)
		population[p2].h = hchild2
		if hchild2 < best_h {
			best_solution = p2
			best_h = hchild2
		}
	}

	return nil
}

func do_mutation() error {
	// fmt.Println("=================================================================")
	parent := rand.Intn(N)

	child, ord := ANNmutation(population[parent])
	if child == nil {
		return nil
	}
	hchildd := calc_tour_lenght(child)

	if !check(child) {
		return errors.New("do_mutation fail")
	}

	if hchildd < population[parent].h {
		count_evol_from_mutation++
		copy(population[parent].chromosome, child)
		copy(population[parent].ord, ord)
		population[parent].h = hchildd
		if hchildd < best_h {
			best_solution = parent
			best_h = hchildd
		}
	}

	return nil
}

func nextGen() error {
	var res error
	b := rand.Float64()
	if b < r {
		res = do_mutation()
	} else {
		res = do_crossover()
	}
	return res

}

func two_opt(solution Solution) int {
	// repeat until no improvement is made
	improve := 0
	for i := 0; i < n-1; i++ {
		for j := i + 1; j < n; j++ {
			new_tour := two_opt_swap(i, j, solution.chromosome)
			if new_tour == nil {
				continue
			}
			new_h := calc_tour_lenght(new_tour)
			if new_h < solution.h {
				improve++
				copy(solution.chromosome, new_tour)
				copy(solution.ord, get_ord(new_tour))
				solution.h = new_h
			}
		}
	}
	return improve
}

func two_opt_swap(i int, k int, tour []int) []int {
	new_tour := make([]int, n)
	// 1. take route[0] to route[i-1] and add them in order to new_route
	copy(new_tour, tour[:i])

	// 2. take route[i] to route[k] and add them in reverse order to new_route
	dec := 0
	for c := i; c < k+1; c++ {
		new_tour[c] = tour[k-dec]
		dec++
	}

	// 3. take route[k+1] to end and add them in order to new_route
	copy(new_tour[k+1:], tour[k+1:])

	if !check(new_tour) {
		fmt.Println(new_tour)
		return nil
	}
	return new_tour
}

func permutations(arr []int) [][]int {
	var helper func([]int, int)
	res := [][]int{}

	helper = func(arr []int, n int) {
		if n == 1 {
			tmp := make([]int, len(arr))
			copy(tmp, arr)
			res = append(res, tmp)
		} else {
			for i := 0; i < n; i++ {
				helper(arr, n-1)
				if n%2 == 1 {
					tmp := arr[i]
					arr[i] = arr[n-1]
					arr[n-1] = tmp
				} else {
					tmp := arr[0]
					arr[0] = arr[n-1]
					arr[n-1] = tmp
				}
			}
		}
	}
	helper(arr, len(arr))
	return res
}

func sort_all_solution() {
	fmt.Println("sort_all_solution, n=", n)
	original := make([]int, n-1)
	for i := 0; i < n-1; i++ {
		original[i] = i + 1
	}
	all_solution := permutations(original)
	n_solution := len(all_solution)
	fmt.Println("n_solution", n_solution)
	min_h := 1e7
	tours := make(ByLength, n_solution)
	for i := 0; i < n_solution; i++ {
		h := calc_tour_lenght2(all_solution[i])
		tours[i] = Edge{i, h}
		if h < min_h {
			min_h = h
		}
	}
	sort.Sort(ByLength(tours))
	for i := 0; i < 5; i++ {
		fmt.Println(i, tours[i].length, "-", all_solution[tours[i].to])
	}
	fmt.Println("min_h", min_h)
}

func ga(init_population InitPopulationFunc, name string) {
	err, best_pos, gpopulation := init_population(N)
	if err != nil {
		fmt.Println("err", err)
		return
	}
	population = gpopulation
	best_h = population[best_pos].h
	best_solution = best_pos
	analyse_edge(true)
	write_edge_count("temp_result/" + name + "_init_edge_count")
	writePopulationH("temp_result/" + name + "_population_after_init")
	fmt.Println("best solution after init", best_h)
	start_time := time.Now()
	h_history := make([]float64, M/100)
	count_evol_from_crossover = 0
	count_evol_from_mutation = 0
	last_count_evol_from_crossover := 0
	last_count_evol_from_mutation := 0
	last_count_evol_from_twoopt := 0

	for idx := 0; idx < M; idx++ {
		if best_h <= best_tour_lenght {
			for i := 0; i < len(h_history); i++ {
				if h_history[i] == 0 {
					h_history[i] = best_h
				}
			}
			break
		}
		err := nextGen()
		if err != nil {
			continue
		}
		// h_history[idx] = best_h
		if idx%100 == 0 {
			h_history[idx/100] = best_h
		}
		if idx%report_period == 0 {
			elapsed := time.Since(start_time)
			count_evol_co := count_evol_from_crossover - last_count_evol_from_crossover
			count_evol_mu := count_evol_from_mutation - last_count_evol_from_mutation
			count_evol_to := count_evol_from_two_opt - last_count_evol_from_twoopt
			fmt.Println("idx", idx, "- count evol", count_evol_co, count_evol_mu, count_evol_to, count_evol_co+count_evol_mu+count_evol_to, " - best_h", best_h, "- elapsed time", elapsed)
			last_count_evol_from_crossover = count_evol_from_crossover
			last_count_evol_from_mutation = count_evol_from_mutation
			last_count_evol_from_twoopt = count_evol_from_two_opt
		}
	}

	elapsed := time.Since(start_time)
	fmt.Println("elapsed time", elapsed)
	count_evol_co := 100. * count_evol_from_crossover / M
	count_evol_mu := 100. * count_evol_from_mutation / M
	fmt.Println("count evol", count_evol_co, "%/", count_evol_mu, "%/", count_evol_mu+count_evol_co)
	fmt.Println("best solution after ga", best_h)
	writePopulationH("temp_result/" + name + "_population_after_ga")
	write_history("temp_result/"+name, h_history)
	write_tour("temp_result/"+name+"_tour", population[best_solution].chromosome)
	// fmt.Println("alpha_value")
	// for i := 0; i < n-1; i++ {
	// 	av := alpha_value[population[best_solution].chromosome[i]][population[best_solution].chromosome[i+1]]
	// 	if av > 0 {
	// 		fmt.Print(av, " ")
	// 	}
	// }
	// av := alpha_value[population[best_solution].chromosome[n-1]][population[best_solution].chromosome[0]]
	// fmt.Println(av)

	analyse_edge(true)
	write_edge_count("temp_result/" + name + "_final_edge_count")
}

func analyse_best_subtour() []int {
	result := make([]int, n)
	for _, chrom := range population {
		tour := chrom.chromosome
		for i := 0; i < n; i++ {
			a := tour[i]
			b := tour[(i+1)%n]
			if best_tour_mask[a] == b || best_tour_mask[b] == a || best_tour_mask[a+n] == b || best_tour_mask[b+n] == a {
				result[i] += 1
			}
		}
	}
	return result
}

func write_edge_count(file_name string) {
	f, _ := os.Create(file_name)
	w := bufio.NewWriter(f)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			w.WriteString(strconv.Itoa(population_edge_mask[i][j]) + "\n")
		}
	}
	w.Flush()
	f.Close()
}

func analyse_edge(reset bool) {
	if reset {
		population_edge_mask = make([][]int, n)
		for i := 0; i < n; i++ {
			population_edge_mask[i] = make([]int, n)
		}
	}
	for _, chrom := range population {
		tour := chrom.chromosome
		for i := 0; i < n-1; i++ {
			population_edge_mask[tour[i]][tour[i+1]] += 1
		}
		population_edge_mask[tour[n-1]][tour[0]] += 1
	}
}

func test_initialization(init_population InitPopulationFunc, name string) {
	h_history := make([]float64, M)
	start_time := time.Now()
	best_subtour_count := make([]int, n)
	avg_best_subtour_count := make([]float64, n)
	population_edge_mask = make([][]int, n)
	for i := 0; i < n; i++ {
		population_edge_mask[i] = make([]int, n)
	}
	for idx := 0; idx < M; idx++ {
		init_population(N)
		population_best_subtour_count := analyse_best_subtour()
		h_history[idx] = best_h
		for i := 0; i < n; i++ {
			best_subtour_count[i] += population_best_subtour_count[i]
			// avg_best_subtour_count[i] = float64(best_subtour_count[i]) / float64(i+1)
		}
		analyse_edge(false)
		if idx%report_period == 0 {
			elapsed := time.Since(start_time)
			fmt.Println("idx", idx, "- elapsed time", elapsed)
		}
	}
	avg_best_subtour_count_n := float64(0)
	for i := 0; i < n; i++ {
		avg_best_subtour_count[i] = float64(best_subtour_count[i]) / float64(M)
		avg_best_subtour_count_n += avg_best_subtour_count[i]
	}
	fmt.Println("avg_best_subtour_count", avg_best_subtour_count_n)
	write_history(name+"_best_subtour_count", avg_best_subtour_count)
	write_history(name+"_test", h_history)
	write_edge_count(name + "_edge_count")
}

func read_population_file() {
	population_file, err := os.Open("nnrawpopulation.txt")
	if err != nil {
		log.Fatal(err)
	}
	defer population_file.Close()
	scanner := bufio.NewScanner(population_file)
	population = make([]Solution, N)
	i := 0
	for scanner.Scan() {
		line := scanner.Text()
		tokens := strings.Fields(line)
		population[i].chromosome = make([]int, len(tokens))
		for j := 0; j < n; j++ {
			gene, _ := strconv.Atoi(tokens[j])
			population[i].chromosome[j] = gene
		}
		// fmt.Println(i, population[i].raw_chromosome)
		population[i].ord = get_ord(population[i].chromosome)
		// fmt.Println(i, population[i].chromosome)
		population[i].h = calc_tour_lenght(population[i].chromosome)

		if population[i].h < best_h {
			// fmt.Println("replace", h, best_solution.h)
			best_solution = i
			best_h = population[i].h
		}

		i++
	}
}

func test_mutation(init_population InitPopulationFunc, name string) {
	err, _, _ := init_population(N)
	if err != nil {
		fmt.Println("err", err)
		return
	}
	parent := rand.Intn(N)
	// fmt.Println("parent", population[parent])
	child, _ := ANNmutation(population[parent])
	write_tour("mutation_test_parent", population[parent].chromosome)
	write_tour("mutation_test_child", child)
}

func gacombine() {
	best_h = 1e10
	population = combined_population
	for i := 0; i < N; i++ {
		if best_h > population[i].h {
			best_h = population[i].h
			best_solution = i
		}
	}
	fmt.Println("population", population[0])
	name := "combined"
	analyse_edge(true)
	write_edge_count(name + "combined_init_edge_count")
	writePopulationH(name + "combined_population_after_init")
	fmt.Println("best solution after init", best_h)
	start_time := time.Now()
	h_history := make([]float64, M)
	count_evol_from_crossover = 0
	count_evol_from_mutation = 0
	for idx := 0; idx < M; idx++ {
		err := nextGen()
		if err != nil {
			break
		}
		h_history[idx] = best_h
		if idx%report_period == 0 {
			elapsed := time.Since(start_time)
			fmt.Println("idx", idx, "- count evol", count_evol_from_crossover, "/", count_evol_from_mutation, "/", count_evol_from_crossover+count_evol_from_mutation, "- best_h", best_h, "- elapsed time", elapsed)
		}
	}
	elapsed := time.Since(start_time)
	fmt.Println("elapsed time", elapsed)
	fmt.Println("count evol", count_evol_from_crossover, "/", count_evol_from_mutation, "/", count_evol_from_crossover+count_evol_from_mutation)
	fmt.Println("best solution after ga", best_h)
	writePopulationH(name + "_population_after_ga")
	write_history(name, h_history)
	write_tour(name+"_tour", population[best_solution].chromosome)

	analyse_edge(true)
	write_edge_count(name + "_final_edge_count")
}

func collect_best(index int, total int) {
	fmt.Println("collect_best", index, "/", total)
	if index == 0 {
		combined_population = make([]Solution, N)
	}
	edges := make(ByLength, N)
	for i := 0; i < N; i++ {
		edges[i] = Edge{i, population[i].h}
	}
	sort.Sort(ByLength(edges))
	start := index * N / total
	end := (index + 1) * N / total
	es := end - start
	fmt.Println("start, end, es", start, end, es)
	for i := 0; i < es; i++ {
		combined_population[i+start].chromosome = make([]int, n)
		combined_population[i+start].ord = make([]int, n)
		copy(combined_population[i+start].chromosome, population[edges[i].to].chromosome)
		copy(combined_population[i+start].ord, population[edges[i].to].ord)
		combined_population[i+start].h = population[edges[i].to].h
		// fmt.Println("i+start", i+start, i)

	}
	fmt.Println("combined_population", combined_population[start].chromosome, combined_population[start].h)
}

func GATS(init_population InitPopulationFunc, name string) {
	for st := 0; st < TS; st++ {
		start_time := time.Now()
		var err error
		err, _, population = init_population(N)
		init_time := time.Since(start_time).Seconds()
		if err != nil {
			fmt.Println("err", err)
			return
		}
		init_h := best_h
		error_rate := (best_h - best_tour_lenght) / best_tour_lenght * 100.
		// fmt.Println("init_h", best_h, "/", best_tour_lenght, " error_rate ", error_rate, best_h-best_tour_lenght)
		poh := make([]float64, N)
		for i, x := range population {
			poh[i] = x.h
		}
		avg_h := stat.Mean(poh, nil)
		avg_conv := (1 - (avg_h-best_tour_lenght)/best_tour_lenght) * 100.
		analyse_edge(true)
		num_edges := 0
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				if population_edge_mask[i][j] > 0 {
					num_edges++
				}
			}
		}
		opt_richness := 0.
		opt_miss_rate := 0.
		phenotype_diversity := float64(num_edges) / float64(n*(n-1)) * 100.
		if opt_tour_avai {
			population_best_subtour_count := analyse_best_subtour()
			num_bst := 0
			num_miss := 0
			for i := 0; i < n; i++ {
				num_bst += population_best_subtour_count[i]
				if population_best_subtour_count[i] == 0 {
					num_miss++
				}
			}
			opt_richness = float64(num_bst) / float64(n*N) * 100.
			opt_miss_rate = float64(num_miss) / float64(n) * 100.

		}
		
		start_time = time.Now()
		for idx := 0; idx < M; idx++ {
			nextGen()
		}
		overall_time := time.Since(start_time).Seconds()+init_time
		f, _ := os.OpenFile("all_results.txt", os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
		w := bufio.NewWriter(f)
		if opt_tour_avai {
			w.WriteString(fmt.Sprintln("Data", file_name, "Alg", name, "init_time", init_time, "overall_time", overall_time,
				"error_rate", error_rate, "avg_conv", avg_conv, "phenotype_diversity", phenotype_diversity,
				"opt_richness", opt_richness, "opt_miss_rate", opt_miss_rate, "init_h", int(init_h), "final_h", int(best_h),
			))
		} else {
			w.WriteString(fmt.Sprintln("Data", file_name, "Alg", name, "init_time", init_time, "overall_time", overall_time,
				"error_rate", error_rate, "avg_conv", avg_conv, "phenotype_diversity", phenotype_diversity,
				"opt_richness", "NA", "opt_miss_rate", "NA", "init_h", int(init_h), "final_h", int(best_h),
			))
		}
		w.Flush()
		f.Close()

	}

}

func gacc() {
	init_populations := []InitPopulationFunc{nn_initilization, random_initilization}
	names := []string{"NN", "RN"}
	crossovers := []DoCrossover{do_crossover, do_crossover_pmx}
	cx_names := []string{"PMX", "GSX1"}
	init_best := float64(0)
	for ip := 0; ip < len(init_populations); ip++ {
		fmt.Println("InitAlg", names[ip])
		for st := 0; st < TS; st++ {
			fmt.Print(st, "-")
			err, _, _ := init_populations[ip](N)
			if err != nil {
				fmt.Println("err", err)
				return
			}
			init_best = best_h
			bkup_population := make([]Solution, N)
			for i := 0; i < N; i++ {
				bkup_population[i].chromosome = make([]int, n)
				copy(bkup_population[i].chromosome, population[i].chromosome)
				bkup_population[i].ord = make([]int, n)
				copy(bkup_population[i].ord, population[i].ord)
				bkup_population[i].h = population[i].h
			}
			for co := 0; co < len(crossovers); co++ {
				var res error
				for idx := 0; idx < M; idx++ {
					if r <= 1 {
						res = crossovers[co]()
						b := rand.Float64()
						if b < r {
							res = do_mutation()
						}
					} else {
						res = do_mutation()
						b := rand.Float64()
						if b < 1/r {
							res = crossovers[co]()
						}
					}
					if err != nil {
						break
					}
				}
				if res != nil {
					fmt.Println(err)
				}
				f, _ := os.OpenFile("compare_results.txt", os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
				w := bufio.NewWriter(f)
				w.WriteString(fmt.Sprintln("Data", file_name, "InitAlg", names[ip], "InitBest", init_best, "cx_names", cx_names[co], "besth", best_h))
				w.Flush()
				f.Close()
				population = bkup_population
				best_h = init_best
			}
		}
		fmt.Println("")

	}
}

func cco() {
	err, _, _ := nn_initilization(N)
	if err != nil {
		fmt.Println("err", err)
		return
	}

	for idx := 0; idx < M; idx++ {
		nextGen()
	}

	sortting := make(ByLength, n)
	for i := 0; i < n; i++ {
		sortting[i] = Edge{i, population[i].h}
	}
	sort.Sort(ByLength(sortting))
	p1 := sortting[0].to
	p2 := sortting[1].to

	// p1 := rand.Intn(N)
	// p2 := rand.Intn(N - 1)
	if p2 >= p1 {
		p2 += 1
	}
	parent1 := population[p1]
	parent2 := population[p2]

	pmx_child1, pmx_child2 := pmx(parent1, parent2)

	p := rand.Intn(n)
	gsx_child1, _ := gsx2(parent1, parent2, p)
	gsx_child2, _ := gsx2(parent2, parent1, p)
	write_tour("temp_result/cco_parent1", parent1.chromosome)
	write_tour("temp_result/cco_parent2", parent2.chromosome)
	write_tour("temp_result/pmx_child1", pmx_child1)
	write_tour("temp_result/pmx_child2", pmx_child2)
	write_tour("temp_result/gsx_child1", gsx_child1)
	write_tour("temp_result/gsx_child2", gsx_child2)
}

func init() {
	arg_n = flag.Int("n", -1, "n")
	arg_file_name = flag.String("file_name", "", "file_name")
	arg_N = flag.Int("N", 100, "N")
	arg_M = flag.Int("M", 100_000, "M")
	arg_TS = flag.Int("TS", 1, "TS")
	arg_r = flag.Float64("r", 0.01, "r")
	arg_T = flag.Float64("T", 50., "T")
	arg_two_opt_rate = flag.Float64("two_opt_rate", 0.01, "two_opt_rate")
	arg_best_tour_length = flag.Float64("BSTTL", 0., "BSTTL")
	arg_init_alg = flag.String("init_alg", "NN", "init_alg")
	arg_report_period = flag.Int("report_period", 10_000, "report_period")
	arg_action = flag.String("action", "GA", "action")
}

func main() {
	flag.Parse()
	n = *arg_n
	file_name = *arg_file_name
	N = *arg_N
	M = *arg_M
	r = *arg_r
	T = *arg_T
	TS = *arg_TS
	two_opt_rate = *arg_two_opt_rate
	init_algs := *arg_init_alg
	report_period = *arg_report_period
	action := *arg_action
	init_algs_list := strings.Split(init_algs, ",")
	fmt.Println("Start, n=", n, "file_name=", file_name)
	fmt.Println("Population =", N)
	fmt.Println("Number of generation =", M)
	fmt.Println("Mutation rate =", r)
	fmt.Println("2-opt rate =", two_opt_rate)
	fmt.Println("T =", T)
	fmt.Println("Init algorithm =", init_algs, init_algs_list)
	fmt.Println("Report period =", report_period)
	fmt.Println("Action =", action)

	rand.Seed(time.Now().UnixNano())
	node_list := read_node_file("../data/" + file_name)
	if node_list == nil {
		fmt.Println("node_list == nil")
		edge_weight = read_edge_weight("../data/" + file_name)
	} else {
		edge_weight = calc_edge_weight(node_list[:])
	}

	calc_nearest()

	for i := 0; i < len(edge_weight); i++ {
		for j := 0; j < len(edge_weight[i]); j++ {
			edge_weight[i][j] = math.Round(edge_weight[i][j])
		}
	}

	best_tour = read_tour_file("../data/" + file_name)
	if best_tour != nil {
		opt_tour_avai = true
		best_tour_lenght = calc_tour_lenght(best_tour[:])
		fmt.Println("Best known tour length", best_tour_lenght)

		best_tour_mask = make([]int, 2*n)
		for i := 0; i < n; i++ {
			best_tour_mask[best_tour[i]] = best_tour[(i+1)%n]
			best_tour_mask[best_tour[i]+n] = best_tour[(i+n-1)%n]
		}
	} else {
		best_tour_lenght = *arg_best_tour_length
		opt_tour_avai = false
	}
	for i, init_alg := range init_algs_list {
		fmt.Println("===============================================================")
		fmt.Println(init_alg)
		init_population := nn_initilization
		switch init_alg {
		case "NN":
			init_population = nn_initilization
		case "MNN":
			init_population = multi_nn_initilization
		case "RNN":
			if !is_nearest_rank {
				calc_nearest_rank()
			}
			init_population = rnn_initilization
		case "MRNN":
			if !is_nearest_rank {
				calc_nearest_rank()
			}
			init_population = multi_rnn_initilization
		case "combined":
			gacombine()
			continue
		default:
			continue
		}
		switch action {
		case "TI":
			test_initialization(init_population, init_alg)
		case "GA":
			ga(init_population, init_alg)
		case "IV":
			sort_all_solution()
		case "GAST":
			GATS(init_population, init_alg)
		case "GAC":
			ga(init_population, init_alg)
			collect_best(i, len(init_algs_list)-2)
		case "GACC":
			gacc()
		case "CCO":
			cco()
		case "TM":
			test_mutation(init_population, init_alg)
		}
	}
}
