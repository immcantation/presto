"""
Multiprocessing functions
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import ctypes
import os
import signal
import sys
import multiprocessing as mp
from collections import OrderedDict
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Defaults import default_coord, default_delimiter, default_out_args
from presto.Annotation import getCoordKey
from presto.IO import getFileType, readSeqFile, countSeqFile, countSeqSets, \
                      getOutputHandle, printLog, printProgress, printWarning, printError

# Constants
TERMINATION_SENTINEL = None
EXCEPTION_SENTINEL = None


class SeqData:
    """
    Class defining sequence data objects for worker processes

    Attributes:
      id : unique identifier
      data : single object or a list of data objects.
      valid : if True data is suitable for processing.
    """
    # Instantiation
    def __init__(self, id, data):
        """
        Initializer

        Arguments:
          id :  unique identifier.
          data : single object or a list of data objects.

        Returns:
          presto.Multiprocessing.SeqData
        """
        self.id = id
        self.data = data
        self.valid = (id is not None and data is not None)

    # Set boolean evaluation to valid value
    def __bool__(self):
        """
        Boolean evaluation

        Returns:
          bool : True if the valid attribute is True
        """
        return self.valid

    # Set length evaluation to number of data records
    def __len__(self):
        """
        Length evaluation

        Returns:
          int : number of objects in the data attribute.
        """
        if isinstance(self.data, SeqRecord) or isinstance(self.data, Seq):
            return 1
        elif self.data is None:
            return 0
        else:
            return len(self.data)


class SeqResult:
    """
    Class defining sequence result objects for collector processes

    Attributes:
      id : unique identifier
      data : single unprocessed object or a list of unprocessed data objects.
      results : single processed object or a list of processed data objects.
      valid : if True processing was successful.
      log : dictionary containing the processing log.
    """
    def __init__(self, id, data):
        """
        Initializer

        Arguments:
          id :  unique identifier.
          data : single unprocessed object or a list of unprocessed data objects.

        Returns:
          presto.Multiprocessing.SeqResult
        """
        self.id = id
        self.data = data
        self.results = None
        self.valid = False
        self.log = OrderedDict([('ID', id)])

    def __bool__(self):
        """
        Boolean evaluation

        Returns:
          bool : True if the valid attribute is True.
        """
        return self.valid

    def __len__(self):
        """
        Length evaluation

        Returns:
          int : number of objects in the results attribute.
        """
        if isinstance(self.results, SeqRecord) or isinstance(self.results, Seq):
            return 1
        elif self.results is None:
            return 0
        else:
            return len(self.results)

    @property
    def data_count(self):
        """
        Data length

        Returns:
          int : number of objects in the data attribute.
        """
        if isinstance(self.data, SeqRecord) or isinstance(self.data, Seq):
            return 1
        elif self.data is None:
            return 0
        else:
            return len(self.data)


def manageProcesses(feed_func, work_func, collect_func,
                    feed_args={}, work_args={}, collect_args={},
                    nproc=None, queue_size=None):
    """
    Manages feeder, worker and collector processes

    Arguments:
      feed_func (function): Data Queue feeder function.
      work_func (function): Worker function.
      collect_func (function): Result Queue collector function.
      feed_args (dict): Dictionary of arguments to pass to feed_func.
      work_args (dict): Dictionary of arguments to pass to work_func.
      collect_args (dict): Dictionary of arguments to pass to collect_func.
      nproc (int): Number of processQueue processes;
                   if None defaults to the number of CPUs
      queue_size (int): Maximum size of the argument queue;
                        if None defaults to 2*nproc

    Returns:
      dict: Dictionary of collector results
    """
    # Define signal handler that raises KeyboardInterrupt
    def _signalHandler(s, f):
        raise SystemExit

    # Define function to terminate child processes
    def _terminate():
        sys.stderr.write('NOTICE> Terminating child processes...  ')
        # Terminate feeders
        feeder.terminate()
        feeder.join()
        # Terminate workers
        for w in workers:
            w.terminate()
            w.join()
        # Terminate collector
        collector.terminate()
        collector.join
        sys.stderr.write('Done.\n')

    # Raise SystemExit upon termination signal
    signal.signal(signal.SIGTERM, _signalHandler)

    # Define number of processes and queue size
    if nproc is None:  nproc = mp.cpu_count()
    if queue_size is None:  queue_size = nproc * 2

    # Define shared child process keep alive flag
    alive = mp.Value(ctypes.c_bool, True)

    # Define shared data queues
    data_queue = mp.Queue(queue_size)
    result_queue = mp.Queue(queue_size)
    # TODO:  find out what's up with this context shenanigans
    ctx = mp.get_context()
    collect_queue = ctx.SimpleQueue()
    # Initiate manager and define shared data objects

    try:
        # Initiate feeder process
        feeder = mp.Process(target=feed_func,
                            args=(alive, data_queue),
                            kwargs=feed_args)
        feeder.start()

        # Initiate worker processes
        workers = []
        for __ in range(nproc):
            w = mp.Process(target=work_func,
                           args=(alive, data_queue, result_queue),
                           kwargs=work_args)
            w.start()
            workers.append(w)

        # Initiate collector process
        collector = mp.Process(target=collect_func,
                               args=(alive, result_queue, collect_queue),
                               kwargs=collect_args)
        collector.start()

        # Wait for feeder to finish and add sentinel objects to data_queue
        feeder.join()
        for __ in range(nproc):  data_queue.put(None)

        # Wait for worker processes to finish and add sentinel to result_queue
        for w in workers:  w.join()
        result_queue.put(None)

        # Wait for collector process to finish and add sentinel to collect_queue
        collector.join()
        collect_queue.put(None)

        # Get collector return values
        collected = collect_queue.get()
    except (KeyboardInterrupt, SystemExit):
        sys.stderr.write('NOTICE> Exit signal received.\n')
        _terminate()
        sys.exit()
    except Exception as e:
        printError('%s.' % e, exit=False)
        _terminate()
        sys.exit()
    except:
        printError('Exiting with unknown exception.', exit=False)
        _terminate()
        sys.exit()
    else:
        if not alive.value:
            printError('Exiting due to child process error.', exit=False)
            _terminate()
            sys.exit()

    return collected


def feedSeqQueue(alive, data_queue, seq_file, index_func=None, index_args={}, ordered=False):
    """
    Feeds the data queue with SeqRecord objects

    Arguments:
      alive : multiprocessing.Value boolean controlling whether processing
              continues; when False function returns
      data_queue : multiprocessing.Queue to hold data for processing
      seq_file : Sequence file to read input from
      index_func : Function to use to define sequence sets
                   if None do not index sets and feed individual records
      index_args : Dictionary of arguments to pass to index_func
      ordered : if True assign sequence indices for deterministic output order

    Returns:
      None
    """
    try:
        # Read input file and index sequence sets if required
        if index_func is None:
            seq_iter = readSeqFile(seq_file)
            if ordered:
                # Assign sequence indices for ordered processing
                data_iter = ((i, s.id, s) for i, s in enumerate(seq_iter))
            else:
                data_iter = ((s.id, s) for s in seq_iter)
        else:
            seq_dict = readSeqFile(seq_file, index=True)
            index_dict = index_func(seq_dict, **index_args)
            if ordered:
                # Assign sequence indices for ordered processing
                data_iter = ((i, k, [seq_dict[j] for j in v]) \
                             for i, (k, v) in enumerate(index_dict.items()))
            else:
                data_iter = ((k, [seq_dict[i] for i in v]) \
                             for k, v in index_dict.items())
    except:
        alive.value = False
        raise

    try:
        # Iterate over data_iter and feed data queue
        while alive.value:
            # Get data from queue
            if data_queue.full():  continue
            else:  data = next(data_iter, None)
            # Exit upon reaching end of iterator
            if data is None:  break

            # Feed queue
            if ordered:
                # Extract sequence index and create SeqData with index
                seq_index, seq_id, seq_data = data
                seq_data_obj = SeqData(seq_id, seq_data)
                seq_data_obj.seq_index = seq_index
                data_queue.put(seq_data_obj)
            else:
                data_queue.put(SeqData(*data))
        else:
            sys.stderr.write('PID %s> Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        raise

    return None


def feedPairQueue(alive, data_queue, seq_file_1, seq_file_2,
                  coord_type=default_coord, delimiter=default_delimiter, ordered=False):
    """
    Feeds the data queue with sequence pairs for processQueue processes

    Arguments:
      alive : a multiprocessing.Value boolean controlling whether processing
              continues; when False function returns
      data_queue : an multiprocessing.Queue to hold data for processing
      seq_file_1 : the name of sequence file 1
      seq_file_2 : the name of sequence file 2
      coord_type : the sequence header format
      delimiter :  a tuple of delimiters for (fields, values, value lists)

    Returns:
      None
    """
    # Function to get coordinate info
    def _key_func(x):
        return getCoordKey(x, coord_type=coord_type, delimiter=delimiter)

    # Generator function to read and check files
    def _read_pairs(seq_file_1, seq_file_2):
        iter_1 = readSeqFile(seq_file_1, index=False)
        iter_2 = readSeqFile(seq_file_2, index=False)
        for seq_1, seq_2 in zip(iter_1, iter_2):
            key_1 = getCoordKey(seq_1.description, coord_type=coord_type,
                                delimiter=delimiter)
            key_2 = getCoordKey(seq_2.description, coord_type=coord_type,
                                delimiter=delimiter)
            if key_1 == key_2:
                yield (key_1, [seq_1, seq_2])
            else:
                raise Exception('Coordinates for sequences %s and %s do not match' \
                                 % (key_1, key_2))

    try:
        # Open and parse input files
        if ordered:
            # For ordered processing, enumerate the pairs to assign indices
            data_iter = enumerate(_read_pairs(seq_file_1, seq_file_2))
        else:
            data_iter = _read_pairs(seq_file_1, seq_file_2)

        # Iterate over data_iter and feed data queue
        while alive.value:
            # Get data from queue
            if data_queue.full():  continue
            else:  data = next(data_iter, None)
            # Exit upon reaching end of iterator
            if data is None:  break

            # Feed queue
            if ordered:
                # Extract sequence index and create SeqData with index
                seq_index, (key, seq_pair) = data
                seq_data_obj = SeqData(key, seq_pair)
                seq_data_obj.seq_index = seq_index
                data_queue.put(seq_data_obj)
            else:
                data_queue.put(SeqData(*data))
        else:
            sys.stderr.write('PID %s> Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        raise

    return None


def processSeqQueue(alive, data_queue, result_queue, process_func, process_args={}):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments:
      alive : multiprocessing.Value boolean controlling whether processing
              continues; when False function returns
      data_queue : multiprocessing.Queue holding data to process
      result_queue : multiprocessing.Queue to hold processed results
      process_func : function to use for processing sequences
      process_args : Dictionary of arguments to pass to process_func

    Returns:
      None
    """
    try:
        # Iterator over data queue until sentinel object reached
        while alive.value:
            # Get data from queue
            if data_queue.empty():  continue
            else:  data = data_queue.get()
            # Exit upon reaching sentinel
            if data is None:  break

            # Perform work
            result = process_func(data, **process_args)
            
            # Propagate sequence index if present for ordered processing
            if hasattr(data, 'seq_index'):
                result.seq_index = data.seq_index

            #import cProfile
            #prof = cProfile.Profile()
            #result = prof.runcall(process_func, data, **process_args)
            #prof.dump_stats('worker-%d.prof' % os.getpid())

            # Feed results to result queue
            result_queue.put(result)
        else:
            sys.stderr.write('PID %s> Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        printError('Error processing sequence with ID: %s.' % data.id, exit=False)
        raise

    return None


def collectSeqQueue(alive, result_queue, collect_queue, seq_file, label,
                    index_field=None, out_file=None, out_args=default_out_args, ordered=False):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments:
      alive : a multiprocessing.Value boolean controlling whether processing
              continues; when False function returns.
      result_queue : Multiprocessing.Queue holding worker results.
      collect_queue : Multiprocessing.Queue to store collector return values.
      seq_file : sample sequence file name.
      label : task label used to tag the output files.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : Common output argument dictionary from parseCommonArgs.
      index_field : Field defining set membership for sequence sets
                    if None data queue contained individual records.
      ordered : if True maintain deterministic output order based on sequence indices.

    Returns:
      None: Adds a dictionary with key value pairs to collect_queue containing
           'log' defining a log object,
           'out_files' defining the output file names
    """
    # Define output format
    out_type = getFileType(seq_file) if out_args['out_type'] is None \
               else out_args['out_type']

    # Wrapper for opening handles and writers
    def _open(x, label=label, out_file=out_file):
        if out_file is not None and x == 'pass':
            handle = open(out_file, 'w')
        else:
            handle = getOutputHandle(seq_file,
                                     out_label='%s-%s' % (label, x),
                                     out_dir=out_args['out_dir'],
                                     out_name=out_args['out_name'],
                                     out_type=out_type)
        return handle

    try:
        # Count records
        if index_field is None:
            result_count = countSeqFile(seq_file)
        else:
            result_count = countSeqSets(seq_file, index_field, out_args['delimiter'])

        # Define log handle
        if out_args['log_file'] is None:
            log_handle = None
        else:
            log_handle = open(out_args['log_file'], 'w')
    except:
        alive.value = False
        raise

    try:
        # Initialize output handles
        pass_handle, fail_handle = None, None
        
        if ordered:
            # For ordered processing, collect all results first, then sort by sequence index
            results_dict = {}
            expected_index = 0
            
            # Collect all results into dictionary keyed by sequence index
            while alive.value:
                # Get result from queue
                if result_queue.empty():  continue
                else:  result = result_queue.get()
                # Exit upon reaching sentinel
                if result is None:  break
                
                # Store result by sequence index
                if hasattr(result, 'seq_index'):
                    results_dict[result.seq_index] = result
                else:
                    # Fallback for results without sequence index
                    results_dict[expected_index] = result
                    expected_index += 1
            
            # Process results in order by sequence index
            start_time = time()
            set_count = seq_count = pass_count = fail_count = 0
            
            for seq_index in sorted(results_dict.keys()):
                result = results_dict[seq_index]
                
                # Print progress for previous iteration
                printProgress(set_count, result_count, 0.05, start_time=start_time)

                # Update counts for current iteration
                set_count += 1
                seq_count += result.data_count

                # Write log
                printLog(result.log, handle=log_handle)

                # Write records
                if result:
                    pass_count += 1
                    try:
                        SeqIO.write(result.results, pass_handle, out_type)
                    except AttributeError:
                        # Open pass file
                        pass_handle = _open('pass')
                        SeqIO.write(result.results, pass_handle, out_type)
                else:
                    fail_count += 1
                    if out_args['failed']:
                        try:
                            SeqIO.write(result.data, fail_handle, out_type)
                        except AttributeError:
                            # Open fail file
                            fail_handle = _open('fail')
                            SeqIO.write(result.data, fail_handle, out_type)
        else:
            # Original unordered processing
            # Iterator over results queue until sentinel object reached
            start_time = time()
            set_count = seq_count = pass_count = fail_count = 0
            while alive.value:
                # Get result from queue
                if result_queue.empty():  continue
                else:  result = result_queue.get()
                # Exit upon reaching sentinel
                if result is None:  break

                # Print progress for previous iteration
                printProgress(set_count, result_count, 0.05, start_time=start_time)

                # Update counts for current iteration
                set_count += 1
                seq_count += result.data_count

                # Write log
                printLog(result.log, handle=log_handle)

                # Write records
                if result:
                    pass_count += 1
                    try:
                        SeqIO.write(result.results, pass_handle, out_type)
                    except AttributeError:
                        # Open pass file
                        pass_handle = _open('pass')
                        SeqIO.write(result.results, pass_handle, out_type)
                else:
                    fail_count += 1
                    if out_args['failed']:
                        try:
                            SeqIO.write(result.data, fail_handle, out_type)
                        except AttributeError:
                            # Open fail file
                            fail_handle = _open('fail')
                            SeqIO.write(result.data, fail_handle, out_type)
            else:
                # Only show error message for unordered processing
                sys.stderr.write('PID %s> Error in sibling process detected. Cleaning up.\n' \
                                 % os.getpid())
                return None

        # Print total counts
        printProgress(set_count, result_count, 0.05, start_time=start_time)

        # Update return values
        log = OrderedDict()
        log['OUTPUT'] = os.path.basename(pass_handle.name) if pass_handle is not None else None
        log['SEQUENCES'] = seq_count
        if index_field is not None:
            log['SETS'] = set_count
        log['PASS'] = pass_count
        log['FAIL'] = fail_count

        # Close file handles and generate return data
        #collect_dict = {'log': log, 'pass': None, 'fail': None}
        collect_dict = {'log': log, 'out_files': []}
        if pass_handle is not None:
            #collect_dict['pass'] = pass_handle.name
            collect_dict['out_files'].append(pass_handle.name)
            pass_handle.close()
        if fail_handle is not None:
            #collect_dict['fail'] = fail_handle.name
            collect_dict['out_files'].append(fail_handle.name)
            fail_handle.close()
        if log_handle is not None:
            log_handle.close()
        collect_queue.put(collect_dict)
    except:
        alive.value = False
        raise

    return None


def collectPairQueue(alive, result_queue, collect_queue, seq_file_1, seq_file_2, label,
                     out_file=None, out_args=default_out_args, ordered=False):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments:
      alive : a multiprocessing.Value boolean controlling whether processing
              continues; when False function returns.
      result_queue : a multiprocessing.Queue holding worker results.
      collect_queue : a multiprocessing.Queue holding collector return values.
      seq_file_1 : the first sequence file name.
      seq_file_2 : the second sequence file name.
      label : task label used to tag the output files.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
      None: adds a dictionary of {log: log object, out_files: output file names} to collect_queue.
    """
    # Define output format
    out_type = getFileType(seq_file_1) if out_args['out_type'] is None \
        else out_args['out_type']

    # Define output names
    if out_args['out_name'] is None:
        out_name_1, out_name_2 = None, None
    else:
        out_name_1 = '%s-1' % out_args['out_name']
        out_name_2 = '%s-2' % out_args['out_name']

    # Wrapper for opening handles and writers
    def _open(x, in_file, out_name, out_file=out_file):
        if out_file is not None and x == 'pass':
            handle = open(out_file, 'w')
        else:
            handle = getOutputHandle(in_file,
                                     out_label='%s-%s' % (label, x),
                                     out_dir=out_args['out_dir'],
                                     out_name=out_name,
                                     out_type=out_type)
        return handle

    try:
        # Count input size
        result_count = countSeqFile(seq_file_1)

        # Define log handle
        if out_args['log_file'] is None:
            log_handle = None
        else:
            log_handle = open(out_args['log_file'], 'w')
    except:
        alive.value = False
        raise

    try:
        # Initialize file handles
        pass_handle, fail_handle_1, fail_handle_2 = None, None, None

        if ordered:
            # For ordered processing, collect all results first, then sort by sequence index
            results_dict = {}
            expected_index = 0
            
            # Collect all results into dictionary keyed by sequence index
            while alive.value:
                # Get result from queue
                if result_queue.empty():
                    continue
                else:
                    result = result_queue.get()
                # Exit upon reaching sentinel
                if result is None:  break
                
                # Store result by sequence index
                if hasattr(result, 'seq_index'):
                    results_dict[result.seq_index] = result
                else:
                    # Fallback for results without sequence index
                    results_dict[expected_index] = result
                    expected_index += 1
            
            # Process results in order by sequence index
            start_time = time()
            iter_count = pass_count = fail_count = 0
            
            for seq_index in sorted(results_dict.keys()):
                result = results_dict[seq_index]
                
                # Print progress for previous iteration
                printProgress(iter_count, result_count, 0.05, start_time=start_time)

                # Update counts for iteration
                iter_count += 1

                # Write log
                printLog(result.log, handle=log_handle)

                # Write assembled sequences
                if result:
                    pass_count += 1
                    try:
                        SeqIO.write(result.results, pass_handle, out_type)
                    except AttributeError:
                        # Open pass file
                        pass_handle = _open('pass', seq_file_1, out_args['out_name'])
                        SeqIO.write(result.results, pass_handle, out_type)
                else:
                    fail_count += 1
                    if out_args['failed']:
                        try:
                            SeqIO.write(result.data[0], fail_handle_1, out_type)
                            SeqIO.write(result.data[1], fail_handle_2, out_type)
                        except AttributeError:
                            # Open fail file
                            fail_handle_1 = _open('fail', seq_file_1, out_name_1)
                            fail_handle_2 = _open('fail', seq_file_2, out_name_2)
                            SeqIO.write(result.data[0], fail_handle_1, out_type)
                            SeqIO.write(result.data[1], fail_handle_2, out_type)
        else:
            # Original unordered processing
            # Iterator over results queue until sentinel object reached
            start_time = time()
            iter_count = pass_count = fail_count = 0
            while alive.value:
                # Get result from queue
                if result_queue.empty():
                    continue
                else:
                    result = result_queue.get()
                # Exit upon reaching sentinel
                if result is None:  break

                # Print progress for previous iteration
                printProgress(iter_count, result_count, 0.05, start_time=start_time)

                # Update counts for iteration
                iter_count += 1

                # Write log
                printLog(result.log, handle=log_handle)

                # Write assembled sequences
                if result:
                    pass_count += 1
                    try:
                        SeqIO.write(result.results, pass_handle, out_type)
                    except AttributeError:
                        # Open pass file
                        pass_handle = _open('pass', seq_file_1, out_args['out_name'])
                        SeqIO.write(result.results, pass_handle, out_type)
                else:
                    fail_count += 1
                    if out_args['failed']:
                        try:
                            SeqIO.write(result.data[0], fail_handle_1, out_type)
                            SeqIO.write(result.data[1], fail_handle_2, out_type)
                        except AttributeError:
                            # Open fail file
                            fail_handle_1 = _open('fail', seq_file_1, out_name_1)
                            fail_handle_2 = _open('fail', seq_file_2, out_name_2)
                            SeqIO.write(result.data[0], fail_handle_1, out_type)
                            SeqIO.write(result.data[1], fail_handle_2, out_type)
            else:
                # Only show error message for unordered processing
                sys.stderr.write('PID %s> Error in sibling process detected. Cleaning up.\n' \
                                 % os.getpid())
                return None

        # Print total counts
        printProgress(iter_count, result_count, 0.05, start_time=start_time)

        # Update return values
        log = OrderedDict()
        log['OUTPUT'] = os.path.basename(pass_handle.name) if pass_handle is not None else None
        log['PAIRS'] = iter_count
        log['PASS'] = pass_count
        log['FAIL'] = fail_count

        # Close file handles and generate return data
        # collect_dict = {'log': log, 'pass': None, 'fail': None}
        collect_dict = {'log': log, 'out_files': []}
        if pass_handle is not None:
            # collect_dict['pass'] = pass_handle.name
            collect_dict['out_files'].append(pass_handle.name)
            pass_handle.close()
        if fail_handle_1 is not None:
            # collect_dict['fail'] = fail_handle.name
            collect_dict['out_files'].append(fail_handle_1.name)
            fail_handle_1.close()
        if fail_handle_2 is not None:
            # collect_dict['fail'] = fail_handle.name
            collect_dict['out_files'].append(fail_handle_2.name)
            fail_handle_2.close()
        if log_handle is not None:
            log_handle.close()
        collect_queue.put(collect_dict)
    except:
        alive.value = False
        raise

    return None
# Convenience functions for ordered processing

def manageProcessesOrdered(feed_func, work_func, collect_func,
                          feed_args={}, work_args={}, collect_args={},
                          nproc=None, queue_size=None):
    """
    Manages feeder, worker and collector processes with ordered output
    
    This is a wrapper around manageProcesses that enables ordered processing
    by adding the 'ordered=True' parameter to feed and collect functions.

    Arguments:
      feed_func (function): Data Queue feeder function.
      work_func (function): Worker function.
      collect_func (function): Result Queue collector function.
      feed_args (dict): Dictionary of arguments to pass to feed_func.
      work_args (dict): Dictionary of arguments to pass to work_func.
      collect_args (dict): Dictionary of arguments to pass to collect_func.
      nproc (int): Number of processQueue processes;
                   if None defaults to the number of CPUs
      queue_size (int): Maximum size of the argument queue;
                        if None defaults to 2*nproc

    Returns:
      dict: Dictionary of collector results
    """
    from presto.Multiprocessing import manageProcesses
    
    # Enable ordered processing in feed and collect functions
    feed_args = dict(feed_args)
    feed_args['ordered'] = True
    
    collect_args = dict(collect_args)
    collect_args['ordered'] = True
    
    return manageProcesses(feed_func, work_func, collect_func,
                          feed_args=feed_args, work_args=work_args, 
                          collect_args=collect_args,
                          nproc=nproc, queue_size=queue_size)


def processSeqFileOrdered(seq_file, process_func, label, 
                         index_func=None, index_args={}, process_args={},
                         out_file=None, out_args=None, 
                         nproc=None, queue_size=None):
    """
    Process sequences from a file with deterministic ordered output
    
    This function provides a high-level interface for processing sequences
    with guaranteed deterministic output order regardless of multiprocessing.

    Arguments:
      seq_file (str): Input sequence file path
      process_func (function): Function to process each sequence/sequence set
      label (str): Label for output files
      index_func (function): Function to group sequences into sets (optional)
      index_args (dict): Arguments for index_func
      process_args (dict): Arguments for process_func
      out_file (str): Output file path (optional, auto-generated if None)
      out_args (dict): Output arguments dictionary
      nproc (int): Number of worker processes (optional)
      queue_size (int): Queue size (optional)

    Returns:
      dict: Dictionary with processing results and output file paths
    """
    from presto.Multiprocessing import feedSeqQueue, processSeqQueue, collectSeqQueue
    from presto.Defaults import default_out_args
    
    if out_args is None:
        out_args = default_out_args
    
    # Determine index field for progress counting
    index_field = None
    if index_func is not None:
        # Try to determine index field from function name or args
        if 'field' in index_args:
            index_field = index_args['field']
    
    # Setup function arguments
    feed_args = {'seq_file': seq_file, 'index_func': index_func, 'index_args': index_args}
    work_args = {'process_func': process_func, 'process_args': process_args}
    collect_args = {'seq_file': seq_file, 'label': label, 'index_field': index_field,
                   'out_file': out_file, 'out_args': out_args}
    
    # Run ordered processing
    return manageProcessesOrdered(feedSeqQueue, processSeqQueue, collectSeqQueue,
                                 feed_args=feed_args, work_args=work_args,
                                 collect_args=collect_args, 
                                 nproc=nproc, queue_size=queue_size)
