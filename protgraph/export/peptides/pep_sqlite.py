import sqlite3
import os
import zlib
from protgraph.export.peptides.abstract_peptide_exporter import \
    APeptideExporter

from protgraph.export.peptides.pep_fasta import PepFasta


class PepSQLite(APeptideExporter):
    """
    A PostGreSQL - Exporter to export PEPTIDES
    into the peptides table

    Those tables will contain all output generated by
    each of the processes. Keep in mind that this table can
    be extremely large, depending on the parmeters set in this tool.

    NOTE: Maybe even exceeding trillion of results for one protein!
    """

    def start_up(self, **kwargs):

        # Traversal parameters:
        self._set_up_taversal(
            kwargs["pep_sqlite_skip_x"],
            kwargs["pep_sqlite_min_pep_length"],
            kwargs["pep_sqlite_miscleavages"],
            kwargs["pep_sqlite_use_igraph"],
            kwargs["pep_sqlite_hops"],
            kwargs["pep_sqlite_batch_size"]
        )
        
        # Create database file
        os.makedirs(kwargs["pep_sqlite_output_dir"], exist_ok=True)
        self.conn = sqlite3.connect(kwargs["pep_sqlite_database"], timeout=10, isolation_level=None)

        self.pepfasta = PepFasta()

        # Create tables if they not exist
        try:
            self._create_tables(**kwargs)
        except Exception as e:
            raise Exception("Could not set up sqlite.", e)

    def _create_tables(self, **kwargs):
        """ Set Pragmas """
        try:
            # Create the large peptides table containing most information
            cur = self.conn.cursor()
            cur.execute("""
            PRAGMA journal_mode=OFF
            """
            )
        except Exception as e:
            print("Warning: Pragma journal_mode was not set, inserting may be slower (Reason: {})".format(str(e)))
        finally:
            self.conn.commit()
            cur.close()

        try:
            # Create the large peptides table containing most information
            cur = self.conn.cursor()
            cur.execute("""
            PRAGMA synchronous=OFF
            """
            )
        except Exception as e:
            print("Warning: pragma synchronous was not set, inserting may be slower (Reason: {})".format(str(e)))
        finally:
            self.conn.commit()
            cur.close()

        try:
            # Create the peptides meta information (can also be extremely large), larger than the peptides tables
            cur = self.conn.cursor()
            cur.execute("""
            CREATE TABLE if not exists peptides_meta (
                peptide TEXT UNIQUE PRIMARY KEY,
                meta TEXT);"""
            )
        except Exception as e:
            print("Warning: Failed creating table 'peptides_meta' (Reason: {})".format(str(e)))
        finally:
            self.conn.commit()
            cur.close()

    def export(self, prot_graph, queue):
        # We continue with the export function
        super().export(prot_graph, queue)

        # and commit everything in the conenction for a protein
        self.conn.commit()

    def export_peptides(self, prot_graph, l_path_nodes, l_path_edges, l_peptide, l_miscleavages, _):
        # Retrieve Meta Infos
        accs = [self.pepfasta._get_accession_or_isoform(prot_graph.vs[nodes[1]]) for nodes in l_path_nodes]
        start_poses = [str(self.pepfasta._get_position_or_isoform_position(prot_graph.vs[nodes[1]])) for nodes in l_path_nodes]
        end_poses = [str(self.pepfasta._get_position_or_isoform_position(prot_graph.vs[nodes[-2]], end=True)) for nodes in l_path_nodes]
        qualifiers =  [",".join(self.pepfasta._get_qualifiers(prot_graph, edges)) for edges in l_path_edges]
        qualifiers = ["," + quali if quali else "" for quali in qualifiers ]
        metas = [
            "".join(
                [
                    acc, "(", str(start_pos), ":", str(end_pos), ",",
                    "mssclvg:", str(misses), quali_entries,  ")"
                ])
            for acc, start_pos, end_pos, quali_entries, misses in zip(accs, start_poses, end_poses, qualifiers, l_miscleavages)
        ]

        # Insert into database
        cur = self.conn.cursor()

        cur.execute('BEGIN TRANSACTION')
        for pep, meta in zip(l_peptide, metas):
            self._retry_query(
                cur, 
                """
                INSERT INTO peptides_meta (peptide, meta)
                VALUES (?, ?)
                ON CONFLICT(peptide) DO UPDATE SET meta = meta || ',' || ? ;
                """,
                [pep, meta, meta]
            )
        cur.execute('COMMIT')
        self.conn.commit()
        cur.close()
        
    def _retry_query(self, cursor, statement, entries):
        # Execute statement. Retry if failed.
        try:
            cursor.execute(statement, entries)
        except Exception:
            self._retry_query(cursor, statement, entries)

    def tear_down(self):
        # Close the connection to postgres
        try:
            self.conn.close()  # Close connection
        except Exception as e:
            print("Connection to Sqlite  could not be closed. (Reason: {})".format(str(e)))
