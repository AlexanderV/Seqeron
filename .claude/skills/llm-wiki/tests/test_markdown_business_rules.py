from __future__ import annotations

import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import wiki_markdown  # noqa: E402


class NavigableLinkContractTests(unittest.TestCase):
    """A link is navigable only when a Markdown reader would see it as prose."""

    def test_escaped_wikilink_is_literal_but_even_backslashes_leave_it_visible(self):
        body = r"\[[escaped]] \\[[visible]] [[ordinary]]"
        self.assertEqual(
            ["visible", "ordinary"],
            wiki_markdown.extract_wikilinks(body),
        )

    def test_blank_target_and_local_anchor_do_not_create_graph_nodes(self):
        body = "[[   ]] [[#installation]] [[guide#installation|Install]]"
        self.assertEqual(["guide"], wiki_markdown.extract_wikilinks(body))

    def test_unmatched_backtick_is_prose_not_an_open_ended_code_span(self):
        body = "An unmatched ` delimiter leaves [[the-link]] visible."
        self.assertEqual(["the-link"], wiki_markdown.extract_wikilinks(body))

    def test_code_span_delimiters_must_have_equal_length(self):
        body = "``code ` [[hidden]] `` [[visible]]"
        self.assertEqual(["visible"], wiki_markdown.extract_wikilinks(body))

    def test_code_span_does_not_cross_heading_or_thematic_break_blocks(self):
        body = "`open\n# [[heading]]\n\n`open\n---\n[[after-break]]"
        self.assertEqual(
            ["heading", "after-break"],
            wiki_markdown.extract_wikilinks(body),
        )

    def test_code_span_does_not_cross_setext_heading_boundary(self):
        body = "`heading\n===\n[[visible]]`"
        self.assertEqual(["visible"], wiki_markdown.extract_wikilinks(body))

    def test_ordered_marker_other_than_one_does_not_interrupt_paragraph(self):
        body = "paragraph\n2.     [[visible-continuation]]"
        self.assertEqual(
            ["visible-continuation"],
            wiki_markdown.extract_wikilinks(body),
        )


class LiteralBlockContractTests(unittest.TestCase):
    """Examples and raw code must never become backlinks or graph edges."""

    def test_wrong_fence_character_and_short_fence_do_not_close_block(self):
        body = (
            "````\n[[hidden-a]]\n~~~\n[[hidden-b]]\n```\n[[hidden-c]]\n"
            "````\n[[visible]]"
        )
        self.assertEqual(["visible"], wiki_markdown.extract_wikilinks(body))

    def test_closing_fence_with_trailing_text_is_not_a_closer(self):
        body = "```\n[[hidden-a]]\n``` trailing\n[[hidden-b]]\n```\n[[visible]]"
        self.assertEqual(["visible"], wiki_markdown.extract_wikilinks(body))

    def test_backtick_in_backtick_fence_info_invalidates_the_opener(self):
        body = "```language`option\n[[visible]]"
        self.assertEqual(["visible"], wiki_markdown.extract_wikilinks(body))

    def test_single_line_html_blocks_and_comments_hide_only_their_own_line(self):
        body = (
            "<pre>[[pre-hidden]]</pre>\n"
            "<!-- [[comment-hidden]] -->\n"
            "[[visible]]"
        )
        self.assertEqual(["visible"], wiki_markdown.extract_wikilinks(body))

    def test_unclosed_html_block_ends_when_its_blockquote_container_ends(self):
        body = "> <script>\n> [[script-hidden]]\n\n[[outside]]"
        self.assertEqual(["outside"], wiki_markdown.extract_wikilinks(body))

    def test_unclosed_fence_ends_when_nested_list_item_ends(self):
        body = (
            "- outer\n"
            "  - inner\n"
            "    ```\n"
            "    [[hidden]]\n"
            "  [[outer-visible]]\n"
        )
        self.assertEqual(["outer-visible"], wiki_markdown.extract_wikilinks(body))

    def test_empty_list_item_still_owns_its_indented_fence(self):
        body = "-\n  ```\n  [[hidden]]\n[[outside]]\n"
        self.assertEqual(["outside"], wiki_markdown.extract_wikilinks(body))

    def test_all_pre_like_html_blocks_hide_wikilinks_until_their_close(self):
        for tag in ("pre", "script", "style", "textarea"):
            with self.subTest(tag=tag):
                body = f"<{tag}>\n[[hidden]]\n</{tag.upper()}>\n[[visible]]"
                self.assertEqual(["visible"], wiki_markdown.extract_wikilinks(body))

    def test_tab_padded_list_continuation_keeps_code_in_the_list(self):
        # The blank line is significant: indented code cannot interrupt the
        # preceding paragraph, even inside a list item.
        body = "- item\n\n\t    [[hidden]]\n\n[[visible]]"
        self.assertEqual(["visible"], wiki_markdown.extract_wikilinks(body))

    def test_crlf_input_preserves_links_outside_code(self):
        body = "[[before]]\r\n~~~\r\n[[hidden]]\r\n~~~\r\n[[after]]\r\n"
        self.assertEqual(
            ["before", "after"],
            wiki_markdown.extract_wikilinks(body),
        )


if __name__ == "__main__":
    unittest.main()
