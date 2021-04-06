import itertools
import multiprocessing

import mysql.connector

from protgraph.export.peptides.abstract_peptide_exporter import \
    APeptideExporter


class PepMySQL(APeptideExporter):
    """
    A MySQL - Exporter to export PEPTIDES
    into the peptides table

    Those tables will contain all output generated by
    each of the processes. Keep in mind that this table can
    be extremely large, depending on the parmeters set in this tool.

    NOTE: Maybe even exceeding trillion of results for one protein!
    """

    @property
    def skip_x(self) -> bool:
        return self.get_skip_x

    @property
    def peptide_min_length(self) -> int:
        return self.get_peptide_min_length

    @property
    def max_miscleavages(self) -> int:
        return self.get_miscleavages

    @property
    def use_igraph(self) -> bool:
        return self.get_use_igraph

    @property
    def peptide_max_length(self) -> int:
        return self.get_peptide_length

    @property
    def batch_size(self) -> int:
        return self.get_batch_size

    def start_up(self, **kwargs):
        # Here we generate a connection to mysql
        # and generate the corresponding tables

        # Connection and other parameters
        self.host = kwargs["pep_mysql_host"]  # Host
        self.port = kwargs["pep_mysql_port"]  # Port
        self.user = kwargs["pep_mysql_user"]  # User
        self.password = kwargs["pep_mysql_password"]  # Password
        self.database = kwargs["pep_mysql_database"]  # Database
        self.no_duplicates = kwargs["pep_mysql_no_duplicates"]

        # Traversal parameters:
        self.get_peptide_length = kwargs["pep_mysql_hops"]  # Number of hops. E.G. 2: s -> h_1 -> h_2 -> e
        self.get_miscleavages = kwargs["pep_mysql_miscleavages"]  # A filter criterion how many miscleavages?
        self.get_peptide_min_length = kwargs["pep_mysql_min_pep_length"]  # Peptide minimum length
        self.get_skip_x = kwargs["pep_mysql_skip_x"]
        self.get_use_igraph = kwargs["pep_mysql_use_igraph"]
        self.get_batch_size = kwargs["pep_mysql_batch_size"]

        # Get Unique generator, since mysql can ONLY! return 1 id at once after bulk inserting entries..
        self.proc_id = multiprocessing.current_process()._identity[0] - 2  # Offset, dur to Main and Reading Thread
        self.id_size = 1000000  # Set sufficiently enough. Fails, if it is called with more than 1 Million processes
        self.num_procs = kwargs["num_of_processes"]
        self.id_gen = self.unique_id_gen()

        # Initialize connection
        try:
            self.conn = mysql.connector.connect(
                user=self.user,
                password=self.password,
                host=self.host,
                database=self.database,
                port=self.port
            )
            # Set a cursor
            self.cursor = self.conn.cursor()
        except Exception as e:
            raise Exception("Could not establish a connection to MySQL (Peptides).", e)

        # Create tables if they not exist
        try:
            self._create_tables(**kwargs)
        except Exception as e:
            raise Exception("Could not create tables in MySQL (Peptides).", e)

    def _create_tables(self, **kwargs):
        """ Create the accessions and peptides tables """
        try:
            # create accessions, so that we only save numbers in the large table!
            cur = self.conn.cursor()
            cur.execute("""
                create table if not exists accessions (
                    id INT NOT NULL AUTO_INCREMENT,
                    accession VARCHAR(15) CHARACTER SET ascii NOT NULL,
                    PRIMARY KEY(id)
                );""")
        except Exception as e:
            print("Error creating accessions table. Continuing... (Reason: {})".format(str(e)))
        finally:
            self.conn.commit()
            cur.close()

        try:
            # Create the large peptides table containing most information
            cur = self.conn.cursor()
            cur.execute("""
            CREATE TABLE  if not exists peptides (
                id {0},
                weight {1} NOT NULL,
                a_count SMALLINT NOT NULL,
                b_count SMALLINT NOT NULL,
                c_count SMALLINT NOT NULL,
                d_count SMALLINT NOT NULL,
                e_count SMALLINT NOT NULL,
                f_count SMALLINT NOT NULL,
                g_count SMALLINT NOT NULL,
                h_count SMALLINT NOT NULL,
                i_count SMALLINT NOT NULL,
                j_count SMALLINT NOT NULL,
                k_count SMALLINT NOT NULL,
                l_count SMALLINT NOT NULL,
                m_count SMALLINT NOT NULL,
                n_count SMALLINT NOT NULL,
                o_count SMALLINT NOT NULL,
                p_count SMALLINT NOT NULL,
                q_count SMALLINT NOT NULL,
                r_count SMALLINT NOT NULL,
                s_count SMALLINT NOT NULL,
                t_count SMALLINT NOT NULL,
                u_count SMALLINT NOT NULL,
                v_count SMALLINT NOT NULL,
                w_count SMALLINT NOT NULL,
                x_count SMALLINT NOT NULL,  -- NOT SKIPPED
                y_count SMALLINT NOT NULL,
                z_count SMALLINT NOT NULL,
                n_terminus VARCHAR(1) CHARACTER SET ascii NOT NULL,
                c_terminus VARCHAR(1) CHARACTER SET ascii NOT NULL,
                PRIMARY KEY (id));""".format(
                    "BINARY(57)" if self.no_duplicates else "BIGINT",
                    "BIGINT" if kwargs["mass_dict_type"] is int else "DOUBLE"
                    )
                )
        except Exception as e:
            print("Error createing peptides table. Continuing... (Reason: {})".format(str(e)))
        finally:
            self.conn.commit()
            cur.close()
            self.peptides_keys = [
                "id", "weight",
                "a_count", "b_count", "c_count", "d_count", "e_count", "f_count", "g_count", "h_count",
                "i_count", "j_count", "k_count", "l_count", "m_count", "n_count", "o_count", "p_count",
                "q_count", "r_count", "s_count", "t_count", "u_count", "v_count", "w_count", "x_count",
                "y_count", "z_count", "n_terminus", "c_terminus"
            ]

        try:
            # Create the peptides meta information (can also be extremely large)
            cur = self.conn.cursor()
            cur.execute("""
            CREATE TABLE  if not exists peptides_meta (
                id BIGINT AUTO_INCREMENT,
                peptides_id {0},
                accession_id INT,
                path MEDIUMTEXT CHARACTER SET ascii NOT NULL,
                miscleavages INT NOT NULL,
                PRIMARY KEY (id)
            );""".format(
                    "BINARY(57)" if self.no_duplicates else "BIGINT",
                )
            )  # References to peptide and accession removed for performance reasons
        except Exception as e:
            print("Error createing peptides_meta table. Continuing... (Reason: {})".format(str(e)))
        finally:
            self.conn.commit()
            cur.close()
            self.peptides_meta_keys = [
                "peptides_id",
                "accession_id",
                "path",
                "miscleavages"
            ]

        # Set insert statement for peptides
        self.statement_accession = "INSERT INTO accessions(accession) VALUES (%s);"
        self.statement_peptides_inner_values = "(" + ",".join(["%s"]*len(self.peptides_keys)) + ")"
        self.statement_peptides_meta_inner_values = "(" + ",".join(["%s"]*len(self.peptides_meta_keys)) + ")"

    def export(self, prot_graph, queue):
        # First insert accession into accession table and retrieve its id:
        # since we only do this per protein!
        accession = prot_graph.vs[0]["accession"]
        self.cursor.execute(
            self.statement_accession,
            (accession,)
        )
        self.accession_id = self.cursor.lastrowid

        # Then we continue with the export function
        super().export(prot_graph, queue)

        # and commit everything in the conenction for a protein
        self.conn.commit()

    def export_peptides(self, prot_graph, l_path_nodes, l_path_edges, l_peptide, l_miscleavages, _):
        # Get the weight
        if "mono_weight" in prot_graph.es[l_path_edges[0][0]].attributes():
            l_weight = [sum(prot_graph.es[x]["mono_weight"]) for x in l_path_edges]
        else:
            l_weight = [-1]*len(l_path_nodes)

        # Generate large list of tuples, which should be added to mysql
        l_peptides_tup = [
            (
                weight,  # Counts of Aminoacids
                peptide.count("A"), peptide.count("B"), peptide.count("C"), peptide.count("D"), peptide.count("E"),
                peptide.count("F"), peptide.count("G"), peptide.count("H"), peptide.count("I"), peptide.count("J"),
                peptide.count("K"), peptide.count("L"), peptide.count("M"), peptide.count("N"), peptide.count("O"),
                peptide.count("P"), peptide.count("Q"), peptide.count("R"), peptide.count("S"), peptide.count("T"),
                peptide.count("U"), peptide.count("V"), peptide.count("W"), peptide.count("X"), peptide.count("Y"),
                peptide.count("Z"),  # N and C Terminus
                peptide[0], peptide[-1]
            )
            for peptide, weight in zip(l_peptide, l_weight)
        ]

        # Insert new entry into database:
        if self.no_duplicates:
            l_peptides_id = self._export_peptide_no_duplicate(l_peptides_tup, l_path_nodes, l_miscleavages)
        else:
            l_peptides_id = self._export_peptide_simple(l_peptides_tup, l_path_nodes, l_miscleavages)

        # Bulk insert meta data information of peptides (ALWAYS)
        l_peptides_meta_tup = [
            (
                peptides_id,
                self.accession_id,
                ",".join(map(str, path_nodes)),
                miscleavages
            )
            for peptides_id, path_nodes, miscleavages in zip(l_peptides_id, l_path_nodes, l_miscleavages)
        ]
        # Bulk insert statement and execute
        stmt = "INSERT INTO peptides_meta (" \
            + ",".join(self.peptides_meta_keys) \
            + ") VALUES " \
            + ",".join([self.statement_peptides_meta_inner_values]*len(l_peptides_id))
        self.cursor.execute(stmt, [y for x in l_peptides_meta_tup for y in x])

    def _export_peptide_simple(self, l_peptides_tup, l_path_nodes, l_miscleavages):
        # Get ids of each simple peptide
        unique_pep_ids = list(itertools.islice(self.id_gen, len(l_peptides_tup)))

        # Bulk Insert them
        stmt = "INSERT INTO peptides (" \
            + ",".join(self.peptides_keys) \
            + ") VALUES " \
            + ",".join([self.statement_peptides_inner_values]*len(l_peptides_tup)) \
            + ";"
        self.cursor.execute(
            stmt,
            [y for a, b in zip(l_peptides_tup, unique_pep_ids) for y in [b] + list(a)]
        )

        # return generated ids
        return unique_pep_ids

    def _export_peptide_no_duplicate(self, l_peptides_tup, l_path_nodes, l_miscleavages):
        # Map peptides to 452 + 4 (multiple of 8) many bits (as id)
        pep_ids = [
            "0000" + "".join(
                # Here we use for the aa-counts 17 bits (each)
                # followed by 5 bits for the n and c terminus (ascii_code - 65)
                # The weight is disregarded, since it is composed by the aa counts
                [format(i, 'b').zfill(17) for i in x[1:-2]] \
                + [format(ord(i) - 65, 'b').zfill(5) for i in x[-2:]]
            )
            for x in l_peptides_tup
        ]

        pep_ids = [int(x, 2).to_bytes(len(x) // 8, byteorder='big') for x in pep_ids]
        # Bulk Insert them
        stmt = "INSERT IGNORE INTO peptides (" \
            + ",".join(self.peptides_keys) \
            + ") VALUES " \
            + ",".join([self.statement_peptides_inner_values]*len(l_peptides_tup)) \
            + ";"
        self._execute_retry(
            stmt,
            [y for a, b in zip(l_peptides_tup, pep_ids) for y in [b] + list(a)]
        )

        # return generated ids
        return pep_ids

    def _execute_retry(self, statement, entries):
        # Execute statement. Retry if failed.
        try:
            self.cursor.execute(statement, entries)
        except Exception:
            self._execute_retry(statement, entries)

    def tear_down(self):
        # Close the connection to mysql
        try:
            self.cursor.close()  # Close cursor
            self.conn.close()  # Close connection
        except Exception as e:
            print("Connection to MySQL  could not be closed. (Reason: {})".format(str(e)))
